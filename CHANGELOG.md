# SNPMeta Changelog
Updates to SNPMeta will be documented here. Versioning follows the form X.Y,
where X is the major version, and Y is the minor version. Major version
increments involve alterations to the way the user interacts with the program,
or fundamentally change how the program behaves. Minor version increments
are small feature additions or bug fixes. Modifications to documentation do not
alter version numbers.


## Old Changes
### 2014-02-19
#### Added
- Counting mechanism, so that SNPMeta moves on if it can't fetch GenBank
  records, for some reason.
- Routine to split the GenBank IDs into smaller lists, so that we
  can avoid a 414 (Request string too long) error.

### 2014-02-14
#### Modified
- Fixed the way SNPMeta looks for SNPs and XML records to annotate
  when the -d and --no-blast options are given.
