Penn Blackbox 
============

The tables in this directory define a simple spherical detector consisting
of 2 R580 PMTs placed at 0 and 90 degrees from a collimated LED source.

`PMTINFO.ratdb` lists the PMT positions and directions.
`blackbox.geo` is the RAT geometry definition.

These and any other experiment-specific database tables in this directory
(e.g. a special noise rate in a DAQ table) are loaded when you put this in
your macro:

    /rat/db/set DETECTOR experiment "blackbox"
    /rat/db/set DETECTOR geo_file "blackbox/blackbox.geo"

