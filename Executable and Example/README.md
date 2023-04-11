**Usage**

Decompress "**modem-exe.zip**" to get"**modem.exe**"(an X64 executable program).
Run the script "**RunForward.Bat**" for forward calculation.

3D Model: **COMMEMI3D-2A.MODmodEM**
Forward Parameter: **COMMEMI3D-2A.modempar**
Model Data: **COMMEMI3D-2A.modemdat**

**Model Data Description**
Add the following 10 field values as input to the original data file of the ModEM:
EPolEx EPolEy EPolHx EPolHy EPolHz
BPolEx BPolEy BPolHx BPolHy BPolHz.
They are added to the file with an initial value of 0 in the same format as the transfer function, and the forward field values can be updated in the file after the program has finished running.