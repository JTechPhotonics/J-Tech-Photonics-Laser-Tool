# J Tech Photonics Laser Tool (Community version)
This Inkscape extension generates gcode for laser cutters and plotting machines from an SVG file.

Version 2.0 just released and there are a lot of changes! If you want you can still access legacy releases (below 2.0) 
 on the [releases page](https://github.com/JTechPhotonics/J-Tech-Photonics-Laser-Tool/releases).
Instructions for older versions can be found on [JTP's website](https://jtechphotonics.com/?page_id=2012).

This extension is essentially a UI wrapped around the [svg_to_gcode](https://github.com/PadLex/SvgToGcode) library. 
So if you want to learn how an Inkscape extension is structured, look no further.
If you're interested in peeking under the hood, check out svg_to_gcode.

## Installation

Download the latest release [here](https://github.com/JTechPhotonics/J-Tech-Photonics-Laser-Tool/releases/latest).
Inkscape versions below 1.0 are not supported. Use legacy releases if you are using Inkscape < 1.0.

Unzip the binaries and copy the laser directly into the Inkscape **user extensions folder**. Inkscape lists the location
of your user extensions folder under **Edit/Inkscape** > **Preferences** > **System**.

Restart Inkscape and you're done.

## Tutorial

### Document Setup
Before using the extension, we need to make sure the document is setup correctly. Open **File** > **Document Properties**.

Set the document's **display units** to `mm` or `in`.
Then set **Scale x**, **Scale y** to `1` and **Viewbox > X**, **Viewbox > Y** to `0`.

<img src="./images/document_setup_properties.png" alt="document_setup_properties.png" width="600" />

Lastly, you can move and rescale your drawing to make it look like it did before. 

### Basic Usage

This extension will parse all svg paths and ignore everything else. 

**Step 1 is to convert all other shapes to paths.** In this case I want to convert the whole drawing to gcode.
So I select everything `ctr+A` and convert the drawing to paths 
**Path** > **Object to Path**.

Open the extension at **Extension** > **Generate Laser Gcode** > **J Tech Community Laser Tool**

Select the **same unit** you used in the **Document Settings**. Then choose an appropriate output directory and 
hit apply.

<img src="./images/important_settings.png" alt="important_settings.png" width="600" />

You'll notice two layers were added to your document:
* `debug reference points` contains the black corners. They 
represent the four corners of your machine's bed. You can use them to eyeball whether the gcode is scaled and placed 
correctly.
* `debug traces` contains the red paths which trace all generated gcode commands.

Note: debug layers are reset everytime you run the extension. So make sure you don't accidentally add any objects to them 
or they will be deleted.

## Contribute

* As a user you can contribute by suggesting features, testing the library and reporting any bugs you encounter in a 
detailed issue.
* As a developer of any skill level you can make pull requests which close issues or introduce useful features. 
Just make sure to create an issue describing what features you want to add before taking the time to implement them.
