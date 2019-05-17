###start with vegetation matrix
###rows=plot numbers X columns=vegetation variables

#creates object in R called veg.table, and attaches it

    veg.table <- read.table("E:\\Kiri Clean\\Veg Analysis\\BecMaster15VegDataSept2018.txt",header=TRUE)
    attach(veg.table)

#produces a summary of veg.table

    summary(veg.table)

#plots all data, saves it

    pdf (file = "E:\\Kiri Clean\\Veg Analysis\\BecMaster15VegDataSept2018.PDF")
    plot(veg.table)
    dev.off()
