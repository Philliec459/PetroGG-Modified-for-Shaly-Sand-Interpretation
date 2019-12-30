# Basic-Well-Log-Interpretation-from-PetroGG-adapted-to-Gulf-Coast-Example

This repository uses Mihai's PetroGG that has been adapted to be used with this Gulf of Mexico Shaley Sand Example.

The Gulf of Mexico data is in GulfCoast.txt file provided. 

PetroGG from Mihai is an excellent foundation for this example. We appreciate all the fine work of Mihai. There have been a few changes made to PetroGG to be more suited to these data. First of all, we have altered the code to use Vshale vs. Vclay. Almost if not all shales have less than 100% clay volume. Therefore, finding endpoint parameters for clay points becomes a challenge since a clay point would become an imaginary point. Instead we use shale point parameters.

Also, we have included two additional Saturation models. Mihai has an excellent Petrophysical repository with various Petrophysical routines including saturations. We started with his Waxman-Smits and Dual-Water saturation models and then made a few changes. For Dual Water we are using the code from George Coates' paper mentioned below. For Waxman-Smits we used the Hill, Shirley and Klein technique to solve for Qv and then used the saturation equation provided by Crain in lieu of an iterative approach.


![Depth_Image](DepthPlot.png)


1. Coates, G.R., Gardner, J.S., and Miller, D.L., 1994, Applying pulse-echo NMR to shaly sand formation evaluation, paper B, 35th Annual SPWLA Logging Symposium Transactions, 22 p.

