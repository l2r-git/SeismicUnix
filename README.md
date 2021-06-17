General concept is to add 3D land and ocean bottom cable geometry handling functionality into Seismic Unix (SU).

This starts with SUTOOLCSV which can input fixed-format text files and parse them into comma-separated-values (CSV) files - with assigned SU key names.
As it parses the files, SUTOOLCSV also examines them and reports errors/issues/problems with the values in the input text files.
In particular, SUTOOLCSV has high-level options which "know" the standard SPS2 and SPS1 files used for 3D and 2D land geometry.
The CSV files output by SUTOOLCSV can be input to SpreadSheets (Microsoft Excell, libreOffice Calc) and repaired/altered/modified if needed due to SU limitations.
SUTOOLCSV can also perform the reverse function (turn these CSV files back into SPS2 or SPS1 fixed-format files).
But SUTOOLCSV is generic and can perform fixed-format to CSV and CSV back to fixed-format for other situations.

Next is SUGEOMCSV which inputs SU seismic traces and updates the SU key values from either fixed-format or CSV files. 
Again, this program has high-level options which "know" the standard SPS2 and SPS1 formats (either as fixed-format files or after conversion to CSV files via SUTOOLCSV).

After that, another program will be written allowing a 3D Grid specification to produce 3D cell numbers (set into the CDP SU key). 
Several other keys will be updated to contain usefull 3D Grid related values. 
Several other computations will occur andbe updated (source-to-receiver OFFSET). 
Further, "nominal" CDP computatuon options will be allowed for 2D survey situations.
