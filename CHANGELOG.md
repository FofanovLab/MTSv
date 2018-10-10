# Change Log

---

## v.1.0.1 - October 4, 2018
### Improved:
**Database and Index Setup**  
- Changed thread prioritization to minimize memory usage.

**Binning and Analysis** 
- Removed special characters from taxa names so output files can be accurately viewed and parsed.
- Sped up and reduced memory required for collecting unaligned queries.

### Fixed:
**Database and Index Setup**
- Fixed memory leak in C++ database building

**Binning and Analysis**
- New signature script fixes issue [#12](https://github.com/FofanovLab/MTSv/issues/12). Taxa with a species group roll up to expected genus and family.

 
### Changed:
**Database and Index Setup**
- Created precursor data file for use in python serialization.

**Binning and Analysis**
- Moved expected value database outside of package to user's $HOME directory so updates do not destroy previous data.
- Unaligned reads are now separated by sample.
- Binning report has more detailed statistics.

### Maintenance:
**Database and Index Setup**
- Stopped using deprecated data structures in C++ code

**Binning and Analysis**
- Added tests for mtst-signature
