# orderly: src

This directory contains a series of source directories for your reports.  Each directory must contain a file `orderly.yml`.  The directory name is important - if you change it it will be treated as a totally separate report.

To create a new report from scratch, use `orderly::orderly_new("name_of_your_report")`. This will be initialised as a new directory `src/name_of_your_report` and will be pre-populated with and `orderly.yml` file. Do not rename this file, just edit it taking careto observe its indentation, blank spaces and punctuation, as they all have meaning. IF in doubt, refer to another `.yml` from an existing (working) report as layout.

PLEASE NOTE: Do not create or edit any orderly reports on master branch!!!! Process should strictly be as follows:

1. Pull (cash) all changes on master
2. Create a new branch, ideally making reference to the orderly report you'll create or edit
3. Work on branch
4. Open a pull request; ensure your branch remains up-to-date with master
5. Test your new/edited orderly report thoroughly!
5. Request a reviewer to approve and merge your branch, once you are sure all works as intended
