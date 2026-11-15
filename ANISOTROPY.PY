
import cv2
import numpy as np
import pandas as pd
import tifffile
import matplotlib.pyplot as plt

from scipy import ndimage as ndi
from skimage.filters import (
    threshold_otsu,
    threshold_triangle,
    threshold_isodata,
    threshold_huang,
    threshold_moments
)
from skimage.segmentation import watershed
from skimage.measure import regionprops, perimeter
from skimage.morphology import disk
from tkinter import Tk, filedialog
import os


# =====================================================
# 1. SELECT ORIGINAL IMAGE
# =====================================================

root = Tk()
root.withdraw()

image_path = filedialog.askopenfilename(
    title="Select Lipid Droplet Image",
    filetypes=[
        ("Image files", "*.tif *.tiff *.png *.jpg *.jpeg"),
        ("All files", "*.*")
    ]
)

if not image_path:
    raise SystemExit("No image selected.")

original = tifffile.imread(image_path) if image_path.lower().endswith(
    (".tif", ".tiff")
) else cv2.imread(image_path, cv2.IMREAD_UNCHANGED)

if original is None:
    raise ValueError("Unable to load image.")

if original.ndim == 3:
    if original.shape[-1] in (3, 4):
        original = cv2.cvtColor(
            original[..., :3], cv2.COLOR_RGB2GRAY
        )
    else:
        raise ValueError("Please select a single-channel image.")

original = np.asarray(original)

# Preserve original intensity for measurement.
measurement_image = original.astype(np.float64)

# Convert to 8-bit only for segmentation.
if original.dtype == np.uint8:
    segmentation = original.copy()
else:
    low, high = np.percentile(original, [0.1, 99.9])

    if high <= low:
        raise ValueError("Image has no intensity variation.")

    segmentation = np.clip(
        (original.astype(float) - low) *
        255.0 / (high - low),
        0, 255
    ).astype(np.uint8)


# =====================================================
# 2. USER SETTINGS
# =====================================================

print("\nLIPID DROPLET RIM INTENSITY ANALYSIS")

methods = ["Otsu", "Triangle", "Default", "Moments", "Huang"]

print("\nThreshold methods:")
for i, method in enumerate(methods):
    print(i + 1, "-", method)

choice = int(input("\nSelect threshold method (1-5): "))
method = methods[choice - 1]

half_width = int(input("Rim half-width (pixels): "))
min_size = int(input("Minimum droplet area (pixels²): "))
max_size = int(input("Maximum droplet area (0 = unlimited): "))
min_circ = float(input("Minimum circularity (0-1): "))

use_watershed = input(
    "Apply watershed? (y/n): "
).lower() == "y"

if half_width < 1 or min_size < 0 or not 0 <= min_circ <= 1:
    raise ValueError("Invalid analysis settings.")

if max_size <= 0:
    max_size = float("inf")


# =====================================================
# 3. THRESHOLDING
# =====================================================

gray = segmentation

threshold_functions = {
    "Otsu": threshold_otsu,
    "Triangle": threshold_triangle,
    "Default": threshold_isodata,
    "Moments": threshold_moments,
    "Huang": threshold_huang
}

threshold_value = threshold_functions[method](gray)

# Assumes bright lipid droplets on a dark background.
binary = gray > threshold_value


# =====================================================
# 4. FILL HOLES
# =====================================================

binary = ndi.binary_fill_holes(binary)


# =====================================================
# 5. WATERSHED SEGMENTATION
# =====================================================

if use_watershed:

    distance = ndi.distance_transform_edt(binary)

    local_max = (
        distance == ndi.maximum_filter(distance, size=3)
    ) & binary

    markers, _ = ndi.label(local_max)

    labels = watershed(
        -distance,
        markers,
        mask=binary
    )

else:

    labels, _ = ndi.label(binary)


# =====================================================
# 6. DETECT INDIVIDUAL DROPLETS
# =====================================================

regions = regionprops(labels)

droplets = []

for region in regions:

    area = region.area

    if area < min_size or area > max_size:
        continue

    droplet_mask = labels == region.label

    perimeter_value = perimeter(droplet_mask)

    if perimeter_value == 0:
        continue

    circularity = (
        4 * np.pi * area / (perimeter_value ** 2)
    )

    if circularity < min_circ:
        continue

    droplets.append({
        "label": region.label,
        "area": area,
        "circularity": circularity,
        "mask": droplet_mask,
        "centroid": region.centroid
    })

n_droplets = len(droplets)

if n_droplets == 0:
    raise SystemExit(
        "No droplets detected. Adjust threshold or size settings."
    )

print("\nDetected droplets:", n_droplets)


# =====================================================
# 7. GENERATE RIM ROI AND MEASURE INTENSITY
# =====================================================

kernel = disk(half_width).astype(np.uint8)

results = []

all_rims = np.zeros(original.shape[:2], dtype=bool)

for i, droplet in enumerate(droplets):

    base_mask = droplet["mask"].astype(np.uint8)

    # Outer boundary
    outer = cv2.dilate(
        base_mask,
        kernel,
        iterations=1
    )

    # Inner boundary
    inner = cv2.erode(
        base_mask,
        kernel,
        iterations=1
    )

    # Rim = outer minus inner
    rim = (outer > 0) & (inner == 0)

    rim_area = np.count_nonzero(rim)

    if rim_area == 0:
        print("Skipping droplet:", i + 1)
        continue

    # Measure on ORIGINAL image
    intensities = measurement_image[rim]

    mean_intensity = np.mean(intensities)
    integrated_density = np.sum(intensities)

    min_intensity = np.min(intensities)
    max_intensity = np.max(intensities)
    std_intensity = np.std(intensities)

    cy, cx = droplet["centroid"]

    results.append({
        "Droplet_ID": i + 1,
        "Droplet_Area": droplet["area"],
        "Circularity": droplet["circularity"],
        "Rim_Area": rim_area,
        "Mean_Rim_Intensity": mean_intensity,
        "Integrated_Density": integrated_density,
        "SD_Rim_Intensity": std_intensity,
        "Min_Intensity": min_intensity,
        "Max_Intensity": max_intensity,
        "Centroid_X": cx,
        "Centroid_Y": cy
    })

    all_rims |= rim


# =====================================================
# 8. SAVE RESULTS
# =====================================================

output_dir = filedialog.askdirectory(
    title="Select Output Folder"
)

if not output_dir:
    raise SystemExit("No output folder selected.")

df = pd.DataFrame(results)

csv_path = os.path.join(
    output_dir,
    "Lipid_Droplet_Rim_Intensity.csv"
)

df.to_csv(csv_path, index=False)

# Save segmentation mask
tifffile.imwrite(
    os.path.join(output_dir, "LD_Segmentation.tif"),
    (labels > 0).astype(np.uint8) * 255
)

# Save individual droplet labels
tifffile.imwrite(
    os.path.join(output_dir, "LD_Labels.tif"),
    labels.astype(np.uint16)
)

# Save combined rim mask
tifffile.imwrite(
    os.path.join(output_dir, "LD_Rim_Mask.tif"),
    all_rims.astype(np.uint8) * 255
)


# =====================================================
# 9. VISUALIZE RESULTS
# =====================================================

plt.figure(figsize=(15, 5))

plt.subplot(1, 3, 1)
plt.imshow(original, cmap="gray")
plt.title("Original Image")
plt.axis("off")

plt.subplot(1, 3, 2)
plt.imshow(labels, cmap="nipy_spectral")
plt.title("Detected Lipid Droplets")
plt.axis("off")

plt.subplot(1, 3, 3)
plt.imshow(original, cmap="gray")
plt.contour(all_rims, levels=[0.5], colors="red")
plt.title("Lipid Droplet Rim ROIs")
plt.axis("off")

plt.tight_layout()
plt.show()


# =====================================================
# 10. FINAL REPORT
# =====================================================

print("\n====================================")
print("ANALYSIS COMPLETED")
print("====================================")
print("Total droplets detected:", n_droplets)
print("Droplets measured:", len(results))
print("Threshold method:", method)
print("Rim half-width:", half_width, "pixels")
print("Approximate rim thickness:", 2 * half_width, "pixels")
print("Results saved at:", csv_path)
print("====================================")

print("\nMean rim intensity:")
print(df["Mean_Rim_Intensity"].describe())
