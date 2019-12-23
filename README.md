# neuracle-eegfile-reader

Neuracle-eegfile-reader runs in Matlab as a plugin of the EEGLAB software. It is used to import Neuracle's EEG data and event files in BDF format.

BDF is a 24 bit version o fthe popular 16 bit EDF format, which was used on previous BioSemi models with 16 bit converters. More details can be found here, <https://www.edfplus.info/>

# Neuracle EEG data files
The Neuracle data format consists of several separate files:

- A binary data.bdf file or together with another data.\*.bdf file which was typically recorded when you checked impedance during recording period, those files containing the voltage values of the EEG  
  
- A binary evt.bdf containing information about events in the data

- A recordInformation.json file containing information about examination such as subject info, exam item and so on.

# Use 

1. Add EEGLAB toolbox folder to the path of Matlab (via Matlab --> File --> Set path)

2. Add neuracle-eegfile-reader to EEGLAB subfolder, plugin folder. 

PS: To import Neuracle's EEG data in EEGLAB, you MUST select one or more data.\*.bdf and the evt.data file at the same time.  

# Licence

```
Copyright (c) 2019 Neuracle, Inc. All Rights Reserved. http://neuracle.cn
```

# Contact
Support Team
- email: support@neuracle.cn
