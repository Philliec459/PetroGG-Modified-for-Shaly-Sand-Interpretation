# Basic-Well-Log-Interpretation-from-PetroGG-adapted-to-Gulf-Coast-Example

This repository started from Mihai's PetroGG that has been adapted to be used with this Gulf of Mexico shaley-sand example.The Gulf of Mexico data is in GulfCoast.txt file provided. 

The foundation for this repository starts with PetroGG from Mihai. PetroGG provides an excellent starting point for our work. We appreciate all the fine work of Mihai on PetroGG. 

We have been a few changes to PetroGG that are well suited to these data. First of all, we have altered the code to employ Vshale vs. Vclay. Almost all shales have less than 100% clay volume. We have personally sampled some of the greasiest, gumbo shales to find that the maximum clay content was only 65% per XRD or FTIR. Therefore, estimating the endpoint parameters for clay points are a challenge since these clay points are imaginary points. Instead we use shale point parameters which are present in most shaley-sands.

Also, we have included two additional Saturation models. Mihai also has an excellent Petrophysical repository with various Petrophysical routines including saturations. We started with his Waxman-Smits and Dual-Water saturation models and then made a few changes. For Dual-Water we are using the code from George Coates (1.) that is primarily used with conventional logs in combinatiomn with NMR data. For Waxman-Smits we use the Hill, Shirley and Klein technique(2.) to solve for Qv from Swb. Finally we then use Waxman-Smits saturation equation provided by Crain in lieu of an iterative approach.

    Qv = Swb/(0.6425/((FluidDensity*Salinity(kppm))**0.5) + 0.22) 

One of the key methods to calibrate the Waxman-Smits results is to refine the m* used in the saturation equation. We do this by plotting the apparent m* vs. Swb to define the lowermost, wet trend in the data. This is the m* trend that defines how m* increases with increasing Swb. m* apparent is calculated using the following equation:

    m* apparent = log10(Rw /(Rt*(1 + Rw*B*Qv))) / log10(PHIT)  

![Mstar_Image](apparent_mstar.png)

Depth Plot:

![Depth_Image](depthPlot.png)


1. Coates, G.R., Gardner, J.S., and Miller, D.L., 1994, Applying pulse-echo NMR to shaly sand formation evaluation, paper B, 35th Annual SPWLA Logging Symposium Transactions, 22 p.

2. Hill, H.J., Shirley, O.J., Klein, G.E.: “Bound Water in Shaley Sands - Its Relation to Qv and Other Formation Properties,” Log Analyst, May-June 1979.


