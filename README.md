# DGating
One-dimensional lattice with gating potential (flashing ratchet)

This is currently the final ratcheting simulation and in my opinion the cleanest.
If you're looking to make your own ratcheting simulation, this is the one to look at right now.
This program was based off of the current version of DLattice, so it's a bit revamped.

If the files have been downloaded, you have to start by giving them permission to run. Type:
  chmod +x ./'Filename'
And then you can start preparing for actually running for data.
You should first look through the GatingStart input data file to set the parameters you'd like to test.
For the first time, A and B can be set in the input file!
Remember that you will only be looking at a particular Force and Gamma in these files, there's not a varying portion for those variables.
Once you're all set, compile the fortran as:
  mpifort -o DGating.exe DGating.f
That's right, the executable is also named better!
Once that's done, run GatingShell.sh and your data should come out after a while.

PQ_Plotter is also the plotting program for the Gating series, with instructions for running that code included inside.
