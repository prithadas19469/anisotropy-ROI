macro "Automatic Lipid Droplet Circles (circle search)" {

    if (nImages == 0)
        exit("Open a fluorescence image first.");

    originalID = getImageID();

    // -------------------------------
    // USER SETTINGS
    // -------------------------------

    Dialog.create("Lipid Droplet Circle Detection");

    Dialog.addNumber("Minimum droplet radius (pixels)", 4);
    Dialog.addNumber("Maximum droplet radius (pixels)", 16);
    Dialog.addNumber("Radius step (pixels, 2 = faster)", 1);
    Dialog.addNumber("Ring thickness (pixels)", 3);
    Dialog.addNumber("Smoothing sigma (pixels, 0 = off)", 0);
    Dialog.addNumber("Detection strength (% of strongest ring, lower = more droplets)", 15);
    Dialog.addNumber("Minimum ring completeness (0-1)", 0.60);
    Dialog.addNumber("Minimum spacing (x sum of radii, lower = allow overlap)", 0.70);

    Dialog.show();

    rMin      = Dialog.getNumber();
    rMax      = Dialog.getNumber();
    rStep     = Dialog.getNumber();
    thickness = Dialog.getNumber();
    sigma     = Dialog.getNumber();
    thrPct    = Dialog.getNumber();
    minComp   = Dialog.getNumber();
    minSpace  = Dialog.getNumber();

    if (rMin < 2 || rMax < rMin || rStep < 1 || thickness <= 0 ||
        sigma < 0 || thrPct <= 0 || thrPct >= 100 ||
        minComp < 0 || minComp > 1 || minSpace <= 0)
        exit("Check your parameters.");

    ringHalf = thickness / 2;
    gap      = 1;      // space between the ring and its inner/outer comparison bands
    flankW   = 3;      // width of the inner/outer comparison bands

    setBatchMode(true);

    // -------------------------------
    // PREPARE A SMOOTHED 32-bit COPY OF THE CURRENT PLANE
    // -------------------------------

    selectImage(originalID);

    if (is("hyperstack")) {
        Stack.getPosition(c, z, t);
        run("Duplicate...", "title=LD_SM duplicate channels=" + c +
            " slices=" + z + " frames=" + t);
    } else {
        run("Duplicate...", "title=LD_SM");
    }

    run("Set Scale...", "distance=0 known=0 unit=pixel");   // work in pixels
    run("32-bit");
    if (sigma > 0)
        run("Gaussian Blur...", "sigma=" + sigma);

    smID = getImageID();
    W = getWidth();
    H = getHeight();

    // -------------------------------
    // RING RESPONSE FOR EVERY RADIUS
    // response = mean(ring band) - mean(inner + outer comparison bands)
    // -------------------------------

    nR = floor((rMax - rMin) / rStep) + 1;
    respID  = newArray(nR);
    radii   = newArray(nR);
    maxID   = 0;

    for (k = 0; k < nR; k++) {

        r = rMin + k * rStep;
        radii[k] = r;
        showStatus("Searching radius " + r + " px ...");
        showProgress(k, nR);

        kern = buildKernel(r, ringHalf, gap, flankW);

        selectImage(smID);
        run("Duplicate...", "title=LDRESP_" + r);
        respID[k] = getImageID();
        run("Convolve...", "text1=[" + kern + "]");

        if (k == 0) {
            run("Duplicate...", "title=LD_MAX");
            maxID = getImageID();
        } else {
            imageCalculator("Max", maxID, respID[k]);
        }
    }

    // -------------------------------
    // FIND CANDIDATE CENTRES
    // -------------------------------

    selectImage(maxID);
    getStatistics(area, mean, min, maxResp);
    thr = maxResp * thrPct / 100;

    nc = 0;
    if (maxResp > 0) {
        run("Find Maxima...", "prominence=" + (thr / 2) + " output=[Point Selection]");
        if (selectionType() != -1) {
            getSelectionCoordinates(cx, cy);
            nc = cx.length;
        }
    }

    // Best radius and score for each candidate centre
    bestScore = newArray(nc);
    bestRad   = newArray(nc);
    for (j = 0; j < nc; j++)
        bestScore[j] = -1e30;

    for (k = 0; k < nR; k++) {
        selectImage(respID[k]);
        for (j = 0; j < nc; j++) {
            v = getPixel(cx[j], cy[j]);
            if (v > bestScore[j]) {
                bestScore[j] = v;
                bestRad[j]   = radii[k];
            }
        }
    }

    // -------------------------------
    // RING COMPLETENESS TEST (on the smoothed copy)
    // counts how many of 24 directions show ring > inner AND ring > outer
    // -------------------------------

    selectImage(smID);

    keep = newArray(nc);
    nKeep = 0;
    off = ringHalf + 2;
    PI2 = 2 * PI;

    for (j = 0; j < nc; j++) {

        if (bestScore[j] < thr) continue;

        x0 = cx[j];
        y0 = cy[j];
        rr = bestRad[j];

        // circle must be fully inside the image
        if (x0 - rr - ringHalf < 0 || y0 - rr - ringHalf < 0 ||
            x0 + rr + ringHalf >= W || y0 + rr + ringHalf >= H)
            continue;

        good = 0;
        for (a = 0; a < 24; a++) {
            ang = a * PI2 / 24;
            ca = cos(ang);
            sa = sin(ang);
            ringV = getPixel(x0 + rr * ca, y0 + rr * sa);
            v2 = getPixel(x0 + (rr - 1) * ca, y0 + (rr - 1) * sa);
            v3 = getPixel(x0 + (rr + 1) * ca, y0 + (rr + 1) * sa);
            if (v2 > ringV) ringV = v2;
            if (v3 > ringV) ringV = v3;
            inV  = getPixel(x0 + (rr - off) * ca, y0 + (rr - off) * sa);
            outV = getPixel(x0 + (rr + off) * ca, y0 + (rr + off) * sa);
            if (ringV > inV && ringV > outV)
                good++;
        }

        if (good / 24 >= minComp) {
            keep[nKeep] = j;
            nKeep++;
        }
    }

    // -------------------------------
    // REMOVE DUPLICATES (strongest first)
    // -------------------------------

    keepScore = newArray(nKeep);
    for (i = 0; i < nKeep; i++)
        keepScore[i] = bestScore[keep[i]];

    order = Array.rankPositions(keepScore);   // ascending
    accepted = newArray(nKeep);
    nAcc = 0;

    for (o = nKeep - 1; o >= 0; o--) {

        j = keep[order[o]];
        ok = true;

        for (q = 0; q < nAcc; q++) {
            jj = accepted[q];
            d = sqrt((cx[j] - cx[jj]) * (cx[j] - cx[jj]) +
                     (cy[j] - cy[jj]) * (cy[j] - cy[jj]));
            if (d < minSpace * (bestRad[j] + bestRad[jj])) {
                ok = false;
                break;
            }
        }

        if (ok) {
            accepted[nAcc] = j;
            nAcc++;
        }
    }

    // -------------------------------
    // CLOSE TEMPORARY IMAGES
    // -------------------------------

    for (k = 0; k < nR; k++) {
        selectImage(respID[k]);
        close();
    }
    selectImage(maxID);
    close();
    selectImage(smID);
    close();

    setBatchMode(false);
    showProgress(1);
    showStatus("");

    selectImage(originalID);

    if (nAcc == 0)
        exit("No droplets detected. Lower 'Detection strength' or 'Minimum ring completeness', " +
             "or check the radius range.");

    // -------------------------------
    // DRAW CIRCLES AS ROIs
    // -------------------------------

    roiManager("Reset");

    for (q = 0; q < nAcc; q++) {
        j  = accepted[q];
        rr = bestRad[j] + ringHalf;
        makeOval(cx[j] - rr, cy[j] - rr, 2 * rr, 2 * rr);
        roiManager("Add");
        roiManager("Select", roiManager("Count") - 1);
        roiManager("Rename", "LD_" + (q + 1));
    }

    // -------------------------------
    // MEASURE ORIGINAL INTENSITY
    // -------------------------------

    selectImage(originalID);

    run("Clear Results");
    run("Set Measurements...",
        "area mean min integrated centroid perimeter redirect=None decimal=3");

    roiManager("Deselect");
    roiManager("Measure");

    // -------------------------------
    // DISPLAY
    // -------------------------------

    roiManager("Deselect");
    roiManager("Show All without labels");
    roiManager("Set Color", "yellow");
    roiManager("Set Line Width", 1);

    print("Lipid droplet analysis completed.");
    print("Detected droplets: " + nAcc);

    showMessage("Analysis completed",
        "Droplets detected: " + nAcc +
        "\nMean intensity is in the Results table.");
}

// ------------------------------------------------------------
// Builds a zero-sum ring-detector kernel as text for Process > Filters > Convolve
//   +1/N on the ring band (radius r +/- ringHalf)
//   -1/M on the inner and outer comparison bands
// ------------------------------------------------------------
function buildKernel(r, ringHalf, gap, flankW) {

    half = Math.ceil(r + ringHalf + gap + flankW);

    innerOut = r - ringHalf - gap;            // inner band: [innerOut - flankW , innerOut)
    outerIn  = r + ringHalf + gap;            // outer band: (outerIn , outerIn + flankW]

    // pass 1: count pixels in each band
    nRing = 0;
    nFlank = 0;
    for (y = -half; y <= half; y++) {
        for (x = -half; x <= half; x++) {
            d = sqrt(x * x + y * y);
            if (abs(d - r) <= ringHalf)
                nRing++;
            else if ((d >= innerOut - flankW && d < innerOut) ||
                     (d > outerIn && d <= outerIn + flankW))
                nFlank++;
        }
    }

    wRing  = 1 / nRing;
    wFlank = -1 / nFlank;

    // pass 2: write the kernel
    s = "";
    for (y = -half; y <= half; y++) {
        row = "";
        for (x = -half; x <= half; x++) {
            d = sqrt(x * x + y * y);
            w = 0;
            if (abs(d - r) <= ringHalf)
                w = wRing;
            else if ((d >= innerOut - flankW && d < innerOut) ||
                     (d > outerIn && d <= outerIn + flankW))
                w = wFlank;
            row = row + d2s(w, 6) + " ";
        }
        s = s + row + "\n";
    }
    return s;
}
