# SNPMeta Changelog
Updates to SNPMeta will be documented here. Versioning follows the form X.Y,
where X is the major version, and Y is the minor version. Major version
increments involve alterations to the way the user interacts with the program,
or fundamentally change how the program behaves. Minor version increments
are small feature additions or bug fixes. Modifications to documentation do not
alter version numbers.

## 2.1 - 2016-10-12
### Modified
- Fixed a bug that caused SNPs with no BLAST hits to either be omitted from the
  final report, or throw a `UnboundLocalError`.

## 2.0 - 2016-02-05
### Modified
- **Transitioned from Python 2 to Python 3.**
- Reorganized script into modules, flow of data through the script should make
  more sense now.
- Fixed a reverse-complement logic error.
- Use temporary files instead of fixed filenames, which can cause problems with
  users' existing files or multiple instances of SNPMeta.

### Added
- Debugging messages to keep track of progress. SNPMeta now prints some
  diagnostic messages related to BLAST searching and record fetching.
- Option to run TBLASTX on local databases.

### Removed
- Verbose flag (`-v`). Output format is now specified by the `--outfmt` option.


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
  when the `-d` and `--no-blast` options are given.
