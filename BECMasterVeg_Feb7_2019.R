###start with vegetation matrix
###rows=plot numbers X columns=vegetation variables

#creates object in R called veg.table, and attaches it

    veg.table <- read.table("I:\\GitHub\\Vegetation_Classification\\BECMasterVeg_Feb7_2019.txt",header=TRUE)
    attach(veg.table)

#produces a summary of veg.table

    summary(veg.table)

#plots all data, saves it

    pdf (file = "I:\\GitHub\\Vegetation_Classification\\BECMasterVeg_Feb7_2019.PDF")
    plot(veg.table)
    dev.off()
