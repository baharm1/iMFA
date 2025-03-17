library(readxl)
library(writexl)

# User-defined inputs ---------------------------------------------------------

# Filenames
datafile = 'input_data_perform_SEnorm/brain_metabolites.xlsx'
factorfile = 'SE_factors.xlsx'
opfilename = 'output_SEnorm/SEnorm_brain_metabolites.xlsx'

# First column  number of MID data
index1 = 4

# Maximum # samples for any one condition
max = 9


# Read and normalize the data -------------------------------------------------

# MID data
metabs = excel_sheets(datafile)
data = lapply(metabs,function(x) data.frame(read_excel(datafile,sheet = x)))
names(data) = metabs

# Normalization factors
factors = read_excel(factorfile)

# Normalize
norm_data = data
for (i in 1:length(metabs))
{
  raw_data = norm_data[[i]]
  ndata = raw_data
  index2 = ncol(raw_data)
  ndata$M0 = NA
  for (j in 1:nrow(raw_data))
  {
    f = factors$factor[factors$Animal_ID == raw_data$Animal_ID[j]]
    f = round(f,2)
    ndata[j,index1:index2] = raw_data[j,index1:index2]*f
    ndata$M0[j] = 100-sum(ndata[j,index1:index2])
  }
  norm_data[[i]] = ndata
}

# Save the data
write_xlsx(norm_data,opfilename)
