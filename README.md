# J Tech Photonics Laser Tool

This is the new official home for the JTP laser cutting Inkscape extension.


## Installation

Download the latest release [here](https://github.com/JTechPhotonics/J-Tech-Photonics-Laser-Tool/releases). Make sure to select the release targeted towards the version of Inkscape you are using. You can find the version of Inkscape you're using under **Help** > **About**.

Unzip the files directly into the Inkscape user extensions folder. Inkscape lists the location of your extensions folder under **Edit** > **Preferences** > **System**.

Restart Inkscape and you're done.


## Documentation

More documentation coming soon. In the meantime you can refer to [JTP's official website](https://jtechphotonics.com/?page_id=2012).

### Custom G-code Header and Footer

Add "header" and "footer" text files without extensions in your destination directory to add custom commands. Don't forget to add a new line at the end of these two files.

If no files are detected the default values are :
- Header : G90 ; Absolute positioning
- Footer : G1 X0 Y0 ; Move to X0 Y0

## For developers

Pull requests are welcome. Just make sure to test your code on both Python >=2.6 and >=3.5 for the main branch, and Python >= 3.6 for the development branch (In the next major release we'll drop support for older python versions in an effort to improve code readability).

