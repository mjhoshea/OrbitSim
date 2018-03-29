 /**
 * Computer Modelling, Exercise 3: 
 *
 * This code is to test the velocity-Verlt algorithm for the Particle3d class.
 *
 * It requires multiple text files to read in the inital states of the Particles and
 * the force constants. 
 *
 * "Initial.input" sets the initial state of both particles and 
 * should contain sequentially for each particle, the first particle is the central one: 
 *
 * one double for the mass
 * six more doubles for position and velocity vectors respectively
 * one string for the particle label
 *
 *
 * @author C. John
 * @author M. O'Shea
 * @version "11/2015"
 */


import java.io.*;

import java.util.Scanner;

public class Particle3dVerlet {


    /*
     * We are using file IO so we need to throw an exception
     */

    /**
     * Main method to test Particle3d class.
     *
     *@param argv[0] name of the output file for position.
     *@param argv[1] name of the output file for the total energy
     *@param argv[2] number of timesteps for time integration
     *@param argv[3] size of each timestep
     *
     */
    
    public static void main (String[] argv) throws IOException{

	//Open output file for position
	String outPos= argv[0];
	PrintWriter positionOutput = new PrintWriter(new FileWriter(outPos));

	//Open output file for energy
	String outEnergy= argv[1];
	PrintWriter energyOutput = new PrintWriter(new FileWriter(outEnergy));
       
	/* Open file containing initial state of central particle
	 * called "Initial.input"
	 *Should contain one double for mass, six more doubles for position and velocity
	 *vectors respectively, one string for particle label
	 */

	
	BufferedReader verletRead = new BufferedReader(new FileReader("Initial.input"));

	//Attach Scanner to file to parse in the information

	Scanner verletScan= new Scanner(verletRead);

	//Create particle from central file

	Particle3d[] mid = new Particle3d[2];
	
	for(int i=0;i<mid.length;i++){
	    mid[i]= new Particle3d();
	}

        mid[0] = Particle3d.readParticle(verletScan);

	mid[1] = Particle3d.readParticle(verletScan);



	
	//Calculate initial force

	Vector3d[] force = new Vector3d[mid.length];
	
	for(int i=0; i<mid.length; i++){
	    force[i]= new Vector3d();
	}

	force = Particle3d.forceCalc(mid); //points from central to orbit
	
	System.out.println(force[0]);
	System.out.println(force[1]);

	//Number of timesteps
	int numstep = Integer.parseInt(argv[2]);

	//Size of timestep
	double dt = Double.parseDouble(argv[3]);

	//Initial time
	double t= 0.0;


	//Test to see if new P.E. method works for arrays

	double potentialE=0;

	potentialE= Particle3d.potentialEnergy(mid);

	
	System.out.println(potentialE);




	//================================================================
	

	// Start of the velocity Verlet algorithm


	//print initial state to position file
	Vector3d posPrint= new Vector3d(mid[1].getPosition());

	positionOutput.printf("%10.5f %10.5f\n", posPrint.getX(), posPrint.getY());

	//print initial total energy to energy file
	double totalE= mid[1].kineticEnergy() + mid[1].potentialEnergy(mid[0]);

	energyOutput.printf("%10.5f %10.7f\n", t, totalE);







	
	for (int i=0;i<numstep;i++){
	
 	    // Update the postion using current velocity and force
	    
	    mid[1].leapPosition(dt, force[1]);
	    
	    // Update forces based on new position

	    Vector3d[] force_new = new Vector3d[mid.length];
	    
	    for(int j= 0; j<mid.length; j++){
		force_new[j]= new Vector3d();
	    }

	    force_new= Particle3d.forceCalc(mid);

	    // Update velocity based on average of current and new force
	    
	    mid[1].leapVelocity(dt,Vector3d.addVector(force[1], force_new[1]).scalarMultiply(0.5));

	    // Set current force to new force

	    force= force_new;
	     
	    // Increase the time
            t = t + dt;

	    //Create vector to assist position printing
	     posPrint= mid[1].getPosition();
 		
	    //Print the current position to file
            positionOutput.printf("%10.5f %10.5f\n", posPrint.getX(), posPrint.getY());

	    //Calculate total energy
	    totalE = mid[1].kineticEnergy() + mid[1].potentialEnergy(mid[0]);

	    //Print the total energy to file
	    energyOutput.printf("%10.5f %10.7f\n", t, totalE);


	}
	// Close the output file
	positionOutput.close();
	


	//Close the energy output file
	energyOutput.close();
    }

}



         
