# CYPRUS

This repository demonstrates Compressive Sensing.

Repository is structured as follows:
```bash
├── src
│   ├── libs
│   │   ├── ...
│   ├── apps
│   │   ├── ...
├── ext
├── docs
├── scripts
```
* src - our code
  * libs - shared code (firmware & software)
  * embedded - focused on embedded systems. ordered by boards, where each board could have the following directories:
    * fw - stores firmware by name of the project (rivers)
    * sw - subdivided based on optional firmwares, where each firmware contains the names of the software apps (trees)
    * hw - board design files
  * apps - pc apps by name of project (flowers)
* ext - external libraries & tool
* docs - documentation, no stricts rules
* scripts - deploy and extract scripts
