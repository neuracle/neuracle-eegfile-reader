# neuracle-eegfile-reader

Neuracle-eegfile-reader runs in Matlab as a plugin of the EEGLAB software. It is used to import Neuracle's EEG data and event files in BDF format.

BDF is a 24 bit version of the popular 16 bit EDF format, which was used on previous BioSemi models with 16 bit converters. More details can be found here, <https://www.edfplus.info/>

# Neuracle EEG data files
There are two kinds of directory structure about recording data which are produced by Neuracle EEG Recorder in the figure below. 

![subject-directory-structure](https://github.com/neuracle/neuracle-eegfile-reader/blob/master/pics/subject-directory-structure.png)

Under *subject_folder* forder, one directory structure contains one or more data.\*.bdf files which are typically created when you check impedance during recording period, an evt.bdf file, and a recorderInfo.json file, the other contains several subfolders where a data.bdf file is included, an evt.bdf file, and a recorderInfo.json file. Sometime other files are also recorded such as audio files and spike.bdf files but not necessary for parsering EEG data. 

- A binary data.\*.bdf file containing the voltage values of the EEG  
  
- A binary evt.bdf file containing information about events in the data

- A recordInformation.json file containing information about examination such as subject info and so on.

# Use 

1. Add EEGLAB toolbox folder to the path of Matlab (via Matlab --> File --> Set path)

2. Add neuracle-eegfile-reader to EEGLAB subfolder, plugin folder.

PS: To import Neuracle's EEG data in EEGLAB, you MUST at the same time select one or more data.\*.bdf files and the evt.data file, or select a folder that contains those files. 

# Licence

```
Copyright (c) 2019 Neuracle, Inc. All Rights Reserved. http://neuracle.cn
```

# Contact
Support Team
- email: support@neuracle.cn
