*************************************************************************

Fluorescence optical diffusion tomography


  Objective:
    - Reconstruct 3-D images of spatially varying fluorescence 
      properties in a turbid media 
      from a frequency-domain optical measurement set
      with an iterative coordinate descent (ICD) optimization 


  System requirements: 
    - Originally developed for Linux 
    - Needs gcc and f77 for source compiling and linking
    - Matlab is used only for processing of the resulting data files. 


  Input/output File formats:
    - Configurations files, such as physical parameter file, 
      and image prior parameter file, are text files. 
    - Image and measurement data is written in a "DAT" file format, 
      which we developed.
      The "MAT" subdirectory contains the matlab M-files 
      for processing this file format. 
    

  References:
  
    - Adam B. Milstein, Jonathan J. Stott, Seungseok Oh, David A. Boas, 
      R. P. Millane, Charles A. Bouman, and Kevin J. Webb, 
      "Fluorescence Optical Diffusion Tomography using Multiple-Frequency Data," 
      Journal of the Optical Society of America A, 
      pp. 1035-1049, vol. 21, no. 6, June 2004. 

    - Adam B. Milstein, Seungseok Oh, Kevin J. Webb, Charles A. Bouman, 
      Quan Zhang, David A. Boas, and R. P. Millane, 
      "Fluorescence Optical Diffusion Tomography,'' 
      Applied Optics, pp 3081-3094, vol. 42, no. 16, June 2003. 

     
  Questions or comments:

    - Prof. Charles Bouman (bouman@purdue.edu)
      Electrical and Computer Engineering, Purdue University

    - Prof. Kevin Webb (webb@purdue.edu)
      Electrical and Computer Engineering, Purdue University

    - Dr. Adam Milstein (amilstei@purdue.edu)
      Lincoln Lab 

	- Justin Patel (patel705@purdue.edu)
      Electrical and Computer Engineering, Purdue University


*************************************************************************

This package consists of 3 subdirectories: 
  1. 'code': inversion source code 
  2. 'example': a test example 
  3. 'getmu': source code used to generate the synthetic absorption 
              and diffusion phantom contained in 'example'
  4. 'getx' : source code used to generate the synthetic fluorescence
              phantom in 'example' directory.
  5. 'MAT': matlab m-files to read/write/process/visualize image or data, 
            which are stored as our specific 'dat' files

*************************************************************************

This package provides 2 scripts to demonstrate how to use it. 
They will show step-by step how to use this package, 
what are the input/output.

  - Type "RUNMEall_with_Matlab_Resolution" in the current directory
    if you have Matlab available in your machine,
    or type "RUNMEall_without_Matlab_Resolution" if you don't. 
    
  (Note)  
  1. What the scripts are doing is 
     1) Compile and link the source code in "code" directory,
        and copy it to ./example/
     2) Generate a simulated measurement data in ./example/
     3) Reconstruct image from the simulated data 
     4) If you are running "RUNMEall_with_Matlab_Resolution", 
        then a sample image will be visualized.
  
*************************************************************************

For programmers who want to modify the sources:

  1. All data structures are defined in "code/structs.h".
  2. All necessary functions are declared in "code/defs.h".

*************************************************************************

Good luck!!!


