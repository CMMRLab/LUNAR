This folder contains the begining works to convert the
Tinker .prm file format to a MSI .frc file format. The
purpose is to get parameters that Tinker offers easily
such that all2lmp.py can provide a larger coverage of FF
types. This work is currently expeirmental and any one
getting this zip may continue this work.

I have not provided the orginal tinker .prm files, so
you will have to go to: https://dasher.wustl.edu/tinker/
and navigate to the Force Field Parameter Sets folder to
download your own versions. This is because Josh Kemppainen
seen a forum where the writer of moltemplate was discussing
that the writer of Tinker would rather you download the files
directly and not send them out in a zip file that is not
orginally tinker.

The generate all2lmp_oplsaa_tinker.frc was converted from
Tinker's oplsaa.prm (not present in this directory, please
download yourself from Tinker), using the read_prm.py and 
the write_oplsaa_frc.py found in this folder.

Josh Kemppainen 11/6/2023 