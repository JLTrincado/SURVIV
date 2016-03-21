## SURVIV: Survival Analysis of mRNA Isoform Variation

Requirements
------------
1. Install Python 2.6.x or Python 2.7.x and corresponding versions of NumPy and
SciPy.
2. Add the Python directory to the $PATH environment variable.

Installation
------------


Usage:
--------------------------------


Example:
--------------------------------
Test SURVIV on sample input:
Run surviv.py as below to test the SURVIV script.

    $ python surviv.py inc.txt surv.txt SURVIV_Result_P.txt

Input: Two input files are required for the SURVIV script. inc.txt: Each row of
this input file contains an alternative splicing event. The 5 columns of this
input file contain the read counts for the two isoforms of the alternative
splicing event. The read counts for different patients are separated by commas
in the column. As an example, for exon skipping events, each row defines a
skipped exon and the columns contain the read counts for inclusion and skipping
isoforms:
- ID: User defined ID for the alternative splicing event.
- IJC: inclusion junction counts, patients are separated by comma.
- SJC: skipping junction counts, patients are separated by comma.
- IncFormLen: length of inclusion form, used for normalization.
- SkipFormLen: length of skipping form, used for normalization.

surv.txt: Each row of this input file contains the survival status for a
patient. Important. The order of the patients in this file should match the
order of patients in inc.txt. The 3 columns of this input file are:
- PatientID: User defined ID for the patient.
- Time: Follow up time.
- Event: The status indicator, 0=alive, 1=dead.

Output: For each alternative splicing event, SURVIV outputs the P-values that
evaluate the associations between alternative splicing and patient survival.
- ID: User defined ID for the alternative splicing event.
- IJC: inclusion junction counts, patients are separated by comma.
- SJC: skipping junction counts, patients are separated by comma.
- IncFormLen: length of inclusion form, used for normalization.
- SkipFormLen: length of skipping form, used for normalization.
- PValue: P-values of the alternative splicing event.


Contacts and bug reports
------------------------
Yi Xing
yxing@ucla.edu

Shihao Shen
shihao@ucla.edu

If you found a bug or mistake in this project, we would like to know about it.
Before you send us the bug report though, please check the following:

1. Are you using the latest version? The bug you found may already have been
   fixed.
2. Check that your input is in the correct format and you have selected the
   correct options.
3. Please reduce your input to the smallest possible size that still produces
   the bug; we will need your input data to reproduce the problem, and the
   smaller you can make it, the easier it will be.

Publication
------------
Shen S, Wu YN, Xing Y. SURVIV: Survival Analysis of mRNA Isoform Variation.
Nature Communications. (In press).

Copyright and License Information
---------------------------------
Copyright (C) 2015 University of California, Los Angeles (UCLA)
Shihao Shen and Yi Xing

Authors: Shihao Shen and Yi Xing

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see http://www.gnu.org/licenses/.
