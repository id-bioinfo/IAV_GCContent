1. run beast xml 'H5G2skygridv2.xml', get .log .trees, run TreeAnotator in BEAST package to get a mcc tree.
2. export rates from Tracer: open .log file in tracer, select statistic to be exported (e.g. location.rates.Austra.Belguim), click 'file' -> 'Export Data Table', export to 'H5G2rates.txt'
3. calculate Bayes factor support for each rate using spread3, input log file, coordinates for each location (coord_dfonly.txt, mostly using centroid of country, can be adjust). https://rega.kuleuven.be/cev/ecv/software/SpreaD3_tutorial#sectionFourTwo. output BF in a text
4. plot in R