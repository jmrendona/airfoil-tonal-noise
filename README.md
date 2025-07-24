# BK-rapid-postproc
Is a collection of functions that allow a rapid post-processing from data obtained by means of the Brüel & Kjær Acoustic adquisition system.

The functions need the have the format given by default by the B&K system. That means that the time should be labeled as 'Ds1-Time', and the pressure from the microphones needes to be such as 'DsX-Signal Y', with X and Y defined values that need to be know in advance by the user of the code. Additionally, for a proper functioning if multiple files want to be analyzed, it is recommended that they have a similar name and only diverge in the script given at the end by the B&K system.  

## Functions included
So far, there are two main functions included, which are those to generate the timetrace and the PSD of the spectra. The images generated are automatically saved and the folders used to store them are also automaticlly created.
For both the inputs are almost the same and read such as:

```
- pr_i = first measurement done with the B&K system (usually 1)
- pr_n = last measurement taken with the same file name. Given by B&K system
- aoa = referes to a variable used to separate the different types of measurements into folders
- chanels = is an array containing the mentioned Y numbers related to the signal channel used in the B&K system
- path = string with the full path where the data is stored. Note that the image will be stored there
- filename = string with the name of the file inclduing the .h5 extention
```

Additionally, to allow a proper comparison of the PSD between all the cases the function fixes the limit of the y-axis with it being an additional user input parameters.

Have fun and good luck with your experiments!
