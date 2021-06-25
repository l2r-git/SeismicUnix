The general idea is to add land and ocean-bottom-cable 3D GEOMETRY handling functionality into Seismic Unix (SU).

This starts with SUTOOLCSV which can input fixed-format text files and parse them into comma-separated-values (CSV) files - with assigned SU key names.
As it parses the files, SUTOOLCSV also examines them and reports errors/issues/problems with the values in the input text files.
In particular, SUTOOLCSV has high-level options which "know" the standard SPS2 and SPS1 files used for 3D and 2D land geometry.
The CSV files output by SUTOOLCSV can be input to SpreadSheets (Microsoft Excell, libreOffice Calc) and repaired/altered/modified if needed to overcome SU limitations.
SUTOOLCSV can also perform the reverse function (turn these CSV files back into SPS2 or SPS1 fixed-format files).
But SUTOOLCSV is generic and can perform fixed-format to CSV and CSV back to fixed-format for other situations and any SU key values.

Next is SUGEOMCSV which inputs SU seismic traces and updates the SU key values from either fixed-format or CSV files. 
Again, this program has high-level options which "know" the standard SPS2 and SPS1 formats (either as fixed-format files or after conversion to CSV files via SUTOOLCSV).
But this program is also generic and can update any SU key values from either fixed-format or CSV files.

After that, another program will be written allowing a 3D Grid specification to produce 3D cell numbers (set into the SU key CDP). 
Several other SU keys will be updated to contain usefull 3D Grid related values. 
Other usefull computations will occur and be updated (source-to-receiver OFFSET). 
Further, a "nominal" 2D CDP computation option will be allowed for 2D survey situations.

The ability to use CSV files for geometry leads in the direction of using SpreadSheets or SQL databases to by-pass the 240 byte header limitations of SU.
Basically you just have to make sure the X, S, and R records of SPS2 format have unique mapping values for their equivalent SpreadSheet/SQL tables produced 
from the CSV files output by SUTOOLCSV - and make sure SUGEOMCSV updates those unique mapping values to known SU trace keys. Which fits easily within 240 bytes.
The mapped SpreadSheet/SQL tables then essentially become external headers to contain information that cannot fit in 240 bytes and also allow processing of
that information without reading the associated SU seismic file.      But this is an extensive concept, and somewhat further in the future.

Also included herein are examples of the 3 SPS2 fixed-format files. And examples of how to run SUTOOLCSV and SUGEOMCSV.
These examples also include problems/issues with the SPS2 files, and possible solutions.

I suggest you review, and run, the examples/tests in the following order:

sutoolcsp_example1.txt
sutoolcsp_errors.txt
sutoolcsp_unrepeat.txt
sutoolcsp_sps2_sps1.txt
sutoolcsp_process.txt

sugeomcsv_create_realistic.txt
sugeomcsv_nicerecord.txt
sugeomcsv_missing.txt
sugeomcsv_scalco_scalel.txt
sugeomcsv_create_ufile.txt
sugeomcsv_statics_and_delim.txt
