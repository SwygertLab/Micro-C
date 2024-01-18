# swygertlab
## Micro Pipeline Github Page

This is a github allocated to providing publically available Micro-C pipelines that will take your Micro-C data and complete the necessary steps in order to analyze Micro-C data. There are 2 versions available. The first version is slightly more simple to use as requires less changes to the overall script and does not require a metadata file, while the second version has a couple additional changes and "requires" a metadata file. However, in exchange for the slight increase in difficulty in version, for which all the necessary steps are still laid out in their respective directories, version 2 also keeps better logs of the steps in the script and has much better organization of the intermediate and final output files.

One thing to note is that Micro-C experiments typically output extremely large datasets. Therefore, these scripts were designed to be run on a system with multiple nodes that users will typically need to understand their own system to modify.
