/**
 *
 * N-body Simulation of the Solar System:
 *
 * This simulation utilizes the velocity-Verlet time-integration method in
 * conjuncture with Newtonian Gravitation to simulate the orbits of 
 * bodies in the solar system.
 *
 *
 * The command line takes five arguments. The first argument is the file name of the file containing the 
 * intial states of every particle. There should be an integer at the start of the file indicating the number
 * of particles in the system. From here, each particle's intial state should be described as follows:
 *
 *=============================================================
 * A string designating the particle label.
 * A double indicating the mass of the particle.
 * Three doubles indicating the initial position of the particle.
 * Three doubles indicating the initial velocity of the particle.
 *============================================================ 
 *
 * Repeat this sequence for all desired particles in the system.
 *
 * IMPORTANT NOTE: The Sun needs to be located at the origin in this input file.
 * The Earth must be the second  particle input.
 * Earth's moon must be the third particle input.
 * 
 *
 * The second command line argument is the file name correlating to the parameters of the system. These
 * should be in the following sequence:
 *
 *===============================================================
 * A double indicating the number of time-steps to evolve over.
 * A double indicating the size of each time-step.
 * A double indicating the value of the gravitational constant, G.
 * =============================================================
 *
 * The third command line argument is the name of the file to write the VMD trajectory info to.
 *
 * The fourth command line argument is the name of the file to write the total energy values to.
 * The fifth command line argument is the name of the file to write the orbit count values and eccentricity values to.
 *
 *
 *
 *@author C. John
 *@author M. O'Shea
 *@version "02/2016"
 */


import java.io.*;
import java.util.Scanner;
import java.text.*;
import java.math.*;

public class OrbitSim {


     /**
     * Main method to run OrbitSim class.
     *
     *@param argv[0] initial state of the particles, with first integer for number of particles in system.
     *@param argv[1] values of the parameters.
     *@param argv[2] name of the output trajectory file for the VMD writing.
     *@param argv[3] name of the output file to monitor energy fluctuations. 
     *@param argv[4] name of the output file to track orbits and eccentricities.
     *
     */

    // We are using file IO so we are throwing an exception

    public static void main (String[] argv) throws IOException{



	/*Generating all necessary file reading and file wriitng capabilities
	 *
	 *
	 *
	 *
	 *
	 */



	//Create reader for initial states of particles
	BufferedReader particleInfo = new BufferedReader(new FileReader(argv[0]));

	//Create scanner from reader for initial states of particles
	Scanner particleScan= new Scanner(particleInfo);

	//Create array of particles from scanner object
	Particle3d[] particles = Particle3d.particleArray(particleScan);

	//Create reader for parameter file
	BufferedReader parameterRead = new BufferedReader(new FileReader(argv[1]));

	//Create scanner from reader for parameter file
	Scanner parameterScan= new Scanner(parameterRead);

	//Create writer for VMD output
        PrintWriter positionOutput = new PrintWriter(new FileWriter(argv[2]));

	//Create writer for energy output
	PrintWriter energyOutput = new PrintWriter(new FileWriter(argv[3]));

	//Create writer for orbit count and eccentricities
	PrintWriter orbitOutput = new PrintWriter(new FileWriter(argv[4]));
   


	
	/* Intialize all parameters, energy, and force array
	 *
	 *
	 *
	 *
	 *
	 */

	
	//Save parameters from scanner object to global variables
	double numStep= parameterScan.nextDouble();
	double sizeStep= parameterScan.nextDouble();
	double G= parameterScan.nextDouble();

	// Set the value of G read from input file
	Particle3d.setG(G);

	//Initial time
	double t=0.0;

	//Create an array of doubles to track aphelions and parahelions
	double[] aphelion= new double[particles.length];
	double[] perihelion= new double[particles.length];

	//Create a double array to track orbits
	double[] orbits = new double[particles.length];

	//Create Vector array of initial positions for orbit tracking
	Vector3d[] posInitial= new Vector3d[particles.length];

	//Create array of position vectors to hold new positions
	Vector3d[] posNext= new Vector3d[particles.length];


	//Create vector for initial separation between Earth and moon
	Vector3d moonInitial= new Vector3d(Vector3d.subVector(Particle3d.seperation(particles[1],particles[2]),particles[1].getPosition()));
	
	//Set-up Vector3d object to hold next seperation between Earth and moon
	Vector3d moonNext= new Vector3d();

	//Create double to track moon orbits around Earth
	double moonOrbits= 0.0;
	
					   


	//initialize  arrays
	for(int i=0;i<particles.length;i++){
	    aphelion[i]= 0.0;
	    perihelion[i]= 100000000;
	    orbits[i]=0;
	    posInitial[i]= particles[i].getPosition();
	    posNext[i]= new Vector3d();
	}


	//Initial energies
	double totalE= Particle3d.potentialEnergy(particles) + Particle3d.kineticEnergy(particles);


	//Calculate Initial Force
	Vector3d[] force = Particle3d.forceCalc(particles);
	

	//Print these initial states to file

	energyOutput.printf("%10.5f %10.5f \n", t, totalE); //something wrong here
	
	Particle3d.toVMD(particles, positionOutput, t);

       

	//Calculate Initial vCoM
	Vector3d vCoM = new Vector3d();

	//Initialize total mass
	double massTotal = 0.0;



	//Adjusting initial velocities of all bodies, such that the CoM
	//remains stationary

	Particle3d.adjustedVelocitys(particles);


	/* Start of loop for time-integration via Velocity-Verlet algorithm
	 *
	 *
	 *
	 *
	 *
	 */

	//Loop for each time step 
	for (int i=0; i<numStep; i++){
	    
	    //test for eccentricities
	    Particle3d.eccentric(perihelion, aphelion, particles);

	    //Leap position of all particles due to current pairwise force  
      	    Particle3d.leapPosition(sizeStep, force, particles);
	    	   

	   
	    //Loop through to save new positions for orbit tracking
	    for(int j=0; j<particles.length; j++){
		posNext[j]= particles[j].getPosition();
	    }
	    
	    //Track orbits
	    Particle3d.orbitTrac(posInitial, posNext, orbits);
	    
	    //Calculate new separation between Earth and moon with corrected
	    moonNext=Vector3d.subVector(Particle3d.seperation(particles[1],particles[2]),particles[1].getPosition());

	    //Track orbits of moon around Earth
	    moonOrbits+= Particle3d.moonTrac(moonInitial, moonNext, moonOrbits);


	    //Update forces based on new positions
	    Vector3d[] forceNew =  Particle3d.forceCalc(particles);
	 
	    //Update velocity based on average of current and new force
	    Particle3d.leapVelocity(sizeStep, force, forceNew, particles);

	    //Set force array to new forces
	    for (int j=0; j<force.length; j++){
	    force[j] = forceNew[j]; 
	    }


	    
	    //Update timestep
	    t += sizeStep;

	    //Print particle positions
	    Particle3d.toVMD(particles, positionOutput,t);

	    //Calculate total energy
	     totalE= Particle3d.potentialEnergy(particles) + Particle3d.kineticEnergy(particles);
	    
	    //Print energy output
	     energyOutput.printf("%10.5f %10.5f \n", t, totalE);

	}
	
    
    
    //Loop to print orbit info
    
    for(int i=0; i<particles.length; i++){
	
	double period= orbits[i]/(numStep*sizeStep);

	    orbitOutput.printf("%s: Perihelion= %10.5f, Aphelion= %10.5f, Orbits= %10.5f, Orbital Period= %10.5f \n\n", particles[i].getLabel(), perihelion[i], aphelion[i], orbits[i], period);
	
	
    }

    //Print number of orbits of moon around Earth to same file
    orbitOutput.printf("Orbits of the moon around Earth: %10.5f \n", moonOrbits);
    
    //Close the output streams
    energyOutput.close();
    positionOutput.close();
    orbitOutput.close();
    
    }
}
