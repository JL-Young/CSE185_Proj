# Origin of Data

Data is taken from CSE185 Lab 6 that was originally taken from the this paper [Functional, metabolic and transcriptional maturation of human pancreatic islets derived from stem cells](https://www.nature.com/articles/s41587-022-01219-z.pdf) and this [GEO website](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5114474). 

# Method
The dataset `GSM5114474_M3_E7` is read in as anndata object. Cells with less than 200 genes expressed and less than 1000 total reads are filtered out. Genes detected in less than 5 cells and genes that have total count of less than 15 is filtered out. 

A random list is generated for the obs and vars. The whole dataset is subsetted by the value of the random number to aroudn 50 obs x 1000 vars. However everytime, the script is ran, the dataset will differ. The anndata is then saved to `h5ad` file format (currently commented out). 

test_data1 has 63 obs x 968 vars.
