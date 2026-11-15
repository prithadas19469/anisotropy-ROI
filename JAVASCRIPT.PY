
macro "Automatic Lipid Droplet Rim ROI" {

    if (nImages == 0)
        exit("Open a fluorescence image first.");

    original = getTitle();

    // -------------------------------
    // USER SETTINGS
    // -------------------------------

    Dialog.create("Lipid Droplet Rim Detection");

    Dialog.addChoice("Threshold method",
        newArray("Otsu", "Triangle", "Default",
                 "Huang", "Moments"), "Otsu");

    Dialog.addNumber("Rim ROI thickness (pixels)", 3);
    Dialog.addNumber("Minimum droplet area (pixels)", 30);
    Dialog.addNumber("Maximum droplet area (0 = unlimited)", 0);
    Dialog.addNumber("Minimum circularity (0-1)", 0.30);
    Dialog.addCheckbox("Apply watershed", false);

    Dialog.show();

    method = Dialog.getChoice();
    thickness = Dialog.getNumber();
    minArea = Dialog.getNumber();
    maxArea = Dialog.getNumber();
    minCirc = Dialog.getNumber();
    watershed = Dialog.getCheckbox();

    if (thickness <= 0 || minArea < 0 ||
        minCirc < 0 || minCirc > 1)
        exit("Check your parameters.");

    if (maxArea <= 0)
        maxArea = 999999999;

    half = thickness / 2;

    // -------------------------------
    // SEGMENTATION
    // -------------------------------

    selectWindow(original);

    run("Duplicate...", "title=LD_MASK");

    selectWindow("LD_MASK");

    run("8-bit");

    setOption("BlackBackground", true);

    setAutoThreshold(method + " dark");

    run("Convert to Mask");

    // Fill holes and connect small gaps
    run("Close-");
    run("Fill Holes");

    if (watershed)
        run("Watershed");

    // -------------------------------
    // DETECT DROPLETS
    // -------------------------------

    roiManager("Reset");

    run("Analyze Particles...",
        "size=" + minArea + "-" + maxArea +
        " circularity=" + minCirc + "-1" +
        " show=Nothing add");

    n = roiManager("Count");

    if (n == 0) {
        close();
        exit("No droplets detected. Adjust threshold.");
    }

    // Store droplet bounding coordinates

    xArray = newArray(n);
    yArray = newArray(n);
    wArray = newArray(n);
    hArray = newArray(n);

    for (i = 0; i < n; i++) {

        roiManager("Select", i);

        getSelectionBounds(x, y, w, h);

        xArray[i] = x;
        yArray[i] = y;
        wArray[i] = w;
        hArray[i] = h;
    }

    // -------------------------------
    // CREATE RIM ROIs
    // -------------------------------

    roiManager("Reset");

    for (i = 0; i < n; i++) {

        x = xArray[i];
        y = yArray[i];
        w = wArray[i];
        h = hArray[i];

        // Skip very small droplets
        if (w <= thickness || h <= thickness)
            continue;

        // Outer ellipse
        makeOval(
            x - half,
            y - half,
            w + thickness,
            h + thickness
        );

        roiManager("Add");

        outerIndex = roiManager("Count") - 1;

        // Inner ellipse
        makeOval(
            x + half,
            y + half,
            w - thickness,
            h - thickness
        );

        roiManager("Add");

        innerIndex = roiManager("Count") - 1;

        // Subtract inner ellipse from outer ellipse
        roiManager("Select",
            newArray(outerIndex, innerIndex));

        roiManager("XOR");

        // Add final ring ROI
        roiManager("Add");

        ringIndex = roiManager("Count") - 1;

        roiManager("Rename",
            "LD_" + (i+1) + "_Rim");

        // Remove temporary outer and inner ROIs
        roiManager("Select",
            newArray(outerIndex, innerIndex));

        roiManager("Delete");
    }

    // -------------------------------
    // REMOVE SEGMENTATION IMAGE
    // -------------------------------

    selectWindow("LD_MASK");
    close();

    // -------------------------------
    // MEASURE ORIGINAL INTENSITY
    // -------------------------------

    selectWindow(original);

    run("Clear Results");

    run("Set Measurements...",
        "area mean integrated centroid perimeter redirect=None decimal=3");

    finalCount = roiManager("Count");

    for (i = 0; i < finalCount; i++) {

        roiManager("Select", i);

        run("Measure");
    }

    // -------------------------------
    // DISPLAY FINAL ROIs
    // -------------------------------

    roiManager("Show All");

    roiManager("Set Color", "yellow");

    roiManager("Set Line Width", 2);

    print("Lipid droplet analysis completed.");
    print("Detected droplets: " + n);
    print("Rim ROI thickness: " + thickness + " pixels");
    print("Measured ROIs: " + finalCount);

    showMessage("Analysis completed",
        "Droplets detected: " + n +
        "\nRim ROIs measured: " + finalCount +
        "\nROI thickness: " + thickness + " pixels");
}
