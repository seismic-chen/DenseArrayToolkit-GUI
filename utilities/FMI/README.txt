FMI readme file

1. What is FMI?
FMI (FISSURES-MATLAB Interface) is an effort to build a bridge between FISSURES
and MATLAB. It allows MATLAB users to take advantage of FISSURES framework to
handle seismic data without knowing the details of FISSURES and Java coding, and
to easily build their own client.

2. Version
This is the version 0.91

3. What are included in this package?
There is a top level directory - FMI, which is the home directory of FMI
package. The files/directories in FMI are:
    @MatEvent           Methods and definition of class MatEvent
    @MatNetwork         Methods and definition of class MatNetwork
    @MatSeismogram      Methods and definition of class MatSeismogram
    apps                GUI applications based on FMI
    FMIdemos            Command-line demos of FMI
    FMItools            Utilities using by FMI
    help                Help files for FMI and GUI
    lib                 Libraries required by FMI
    m_map               MATLAB Mapping toolbox
    matTaup             MATLAB version of TauP toolbox
    README.txt          FMI readme
    install.txt         FMI installation instruction

4. System Requirements
FMI can be run on any system with MATLAB version 6.5 (R13) or high installed.
The Java Virtual Machine (JVM) must be enabled. A system with 500MHz CPU or
above, and 256 MB memory or above is recommended.

5. Dependencies
1) FISSURES framework. The following Java libraries are required
	fissuresIDL       1.0
	fissuresImpl      1.1.4
	fissuresUitl      1.0.6  (namingService only)
	SeedCodec         1.0Beta2
	JEvalRespClasses  1.2.01
	matTaupClasses    1.0
	taup              1.1.4
	junit             3.8.1
	junit-addon       1.4
	OB                4.1.0
	OBNaming          4.1.0

Please note that all required Java classes are packed into a single archive
file FMI.jar. You only need to add FMI.jar into classpath.txt in MATLAB.

2) m_map toolbox. The standard m_map toolbox only includes a very coarse line
data file. The high-resolution GSHHS coastline data files are optional. You can
download them from FMI website or from University of Hawaii website 
(http://www.soest.hawaii.edu/wessel/gshhs/gshhs.html).

3) matTaup toolbox. A set MATLAB tools to utilize functions in TauP toolbox. 
Running matTaup needs Java Taup toolkit and matTaupClasses, which are already
included in FMI.jar


6. Installation and Usage
Please follow the instructions in install.txt. After installation, read
FMI.doc in FMIHOME/help directory for details of FMI programming. Several demos
and two example applications based on FMI are also provided.

7. License
This an open-source software. The FMI project is developed and released under 
the terms of Gnu Public License (http://www.gnu.org/copyleft/gpl.html).

8. Questions and bugs?
For project management, contact Ken Creager (kcc@ess.washington.edu)
For FMI bugs, contact Qin Li (qinli@u.washington.edu)
For EventFinder and SeisFinder, contact Ronnie Ning (ronnie@ess.washington.edu) 

You can also send any questions and comments to fmihelp@ess.washington.edu

