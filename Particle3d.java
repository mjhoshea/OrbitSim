/**
 * Particle3D class.
 *
 * This class utilises the properties of Vector3d class
 * in order to create arrays of particles. Contained within the class are numerous
 * methods to carry out calculations. These possess the functionality to set up
 * an array of particles and format there info in VMD format,
 * time-integrate positions and velocity's, calculate the potential and kinetic energies 
 *of a system of particles, calculate pairwise Gravitational forces between particles, calculate
 * the eccentricity and finally track the total orbits in a simulation. 
 *
 *
 * @auther M. OSHEA
 * @auther C. JOHN
 * @version "11/2014"
 *
 */

import java.util.Scanner;
import java.lang.Math;
import java.io.*;



public class Particle3d {
    
    /* ******************************************
     * Properties
     ********************************************/
    
    private double mass;
    private Vector3d position;
    private Vector3d velocity; 
    private String label;  
    private static double G;
    
    // Setters and Getters

    /** Get the position of a particle.
     *
     * @return a Vector3d representing the position.
     */

    public Vector3d getPosition() { return position; }

    /** Get the velocity of a particle.
     *
     * @return a Vector3d representing the velocity.
     */

    public Vector3d getVelocity() { return velocity; }
    
    /** Get the mass of a particle.
     *
     * @return a double representing the mass.
     */
    
    public double getMass() { return mass; }
    
    /** Get the label of a particle.
     *
     *@return a string representing the label.
     */

    public String getLabel() {return label;}
       
    /** Get the value of G
     *
     *
     *@return a double that is G
     */

    public static double getG(){return G;}


    /** Set the value of G
     *
     *@param val a double to set the value of G to.
     */

    public static void setG(double val){ G=val;} 



    /** Set the position of a particle.
     *
     * @param p a Vector3d representing the position.
     */
    
    public void setPosition(Vector3d p) { this.position = p;}

    /** Set the velocity of a particle.
     *
     * @param v a Vector3d representing the velocity.
     */
    
    public void setVelocity(Vector3d v) { this.velocity = v;}
    
    /** Set the mass of a particle.
     *
     * @param m a double representing the mass.
     */
    
    public void setMass(double m) { this.mass = m; } 


    /** Set the label of a particle.
     *
     *@param l a string representing the label.
     */

    public void setLabel(String l){this.label= l;}


    
    /* ******************************************
     * Constructors
     ********************************************/
    
    /** Default constructor. Sets all properties to a value
     * of zero.
     */
    
    public Particle3d() {
	this.setMass(0); 
	this.setPosition(new Vector3d());
	this.setVelocity(new Vector3d());
	this.setLabel("na");
    }
    /** Explicit constructor. Constructs a new Particle1D with
     * explicitly given position, velocity, and mass.
     *
     * @param m a double that defines the mass.
     * @param p a double that defines the position.
     * @param v a double that defines the velocity.
     * @param l a String that names the particle  
     */
    
    public Particle3d(double m, Vector3d p,  Vector3d v, String l) {
	this.setMass(m);
	this.setPosition(p);
	this.setVelocity(v);
	this.setLabel(l);
    }
    
    /* ******************************************
     * toString Method
     ********************************************/
    
    
    /** Returns a String representation of Particle3d.
     * Used to print a Particle3d instance using the "%s"
     * format identifier.
     */
    public String toString() {
	return  this.getLabel()  + " " + this.getPosition()+" " ;
    }
    
    /* ******************************************
     * Instance Methods
     ********************************************/
    
    
    /** Time integration support: evolve the position
     * according to dx = v * dt.
     *
     * @param dt a double that is the timestep.
     */
    public void leapPosition(double dt) {
	
        position = Vector3d.addVector(position, velocity.scalarMultiply(dt));
    }



    /** Time integration support: evolve the position
     * according to dx = x + v * dt + 0.5 * a * dt^2.
     *
     * @param dt a double that is the timestep.
     * @param force an array of vectors that stores the current force.
     * @param particles that is an array of Particle3d objects.
     */
    
    public static void leapPosition(double dt, Vector3d[] force, Particle3d[] particles){
	
	//Initialise acceleration as a Vector3d and set values to zero  
	Vector3d acceleration = new Vector3d();
	
	//Loop to cycle through all particles in array
	
	for (int i = 0 ; i < particles.length ; i ++){
	    
	    //Calculates the acceleration for each particle in the array and uses 
	    //this to calculate the updated position
	    
	    acceleration = force[i].scalarDivide(particles[i].getMass());
	    
	    Vector3d  posI = new Vector3d(particles[i].getPosition());
	    Vector3d  velI= new Vector3d(particles[i].getVelocity());
	    
	    particles[i].setPosition( Vector3d.addVector(Vector3d.addVector(posI, velI.scalarMultiply(dt)), acceleration.scalarMultiply(0.5*dt*dt)));
   
    }
    
   }

    
   /** Time integration support: evolve the velocity
    * according to dv = (f1+f2)/2m * dt.
    *
    * @param dt a double that is the timestep.
    * @param force an array of vectors that stores the current force.
    * @param particles that is an array of Particle3d objects.
    * @param forceNew an array of Vector3d objects to store updated force
    */
   
    public static void leapVelocity(double  dt, Vector3d[] force, Vector3d[] forceNew, Particle3d[] particles) {
	
	Vector3d updateV = new Vector3d();
      
	//Loop to cycle through all particles in array

    	for (int i = 0; i < particles.length ; i ++){

	    //initialize mass variable
	    double mass= particles[i].getMass();

	    
	    updateV = Vector3d.addVector(force[i].scalarDivide(2*mass), forceNew[i].scalarDivide(2*mass));
    	 particles[i].setVelocity( Vector3d.addVector(particles[i].getVelocity(), updateV.scalarMultiply(dt)));
        

	}
    }



    /** Returns the kinetic energy of a Particle3d,
     * calculated as 1/2*m*v^2. Where v, is the magnitude of the
     *velocity vector . 
     *
     * @return a double that is the kinetic energy.
     */
    
    public double kineticEnergy() { return 0.5*mass*velocity.mag()*velocity.mag();}


    
    /* ******************************************
     * Static Methods
     ********************************************/


    /** Method to calculate total kinetic energy of a system of particles
     * stored within an array. Calculations performed utilizing the 
     *kineticEnergy instance method described above.
     *
     *
     *
     *@param particles the array containing all particles of the system.
     *
     *@return double that contains the sum of individual kinetic energy 
     * contributions from each particle.
     */

public static double kineticEnergy(Particle3d[] particles){
	
    double sum = 0; 

    for(int i=0;i<particles.length;i++){

	sum += particles[i].kineticEnergy();
    }

    return  sum;
}



    
    /** Method to calculate Gravitational Potential energy of two particles
     *
     *@param particles an array of Particle3d objects containing all particles in the system.
     *
     *
     *@return double that is the potential energy of the system.
     */


    public static double potentialEnergy(Particle3d[] particles){

	double potential= 0;

	
	//Loop to cycle through particles for energy calculation
	for(int i=0;i<particles.length;i++){
	   
	    //Loop to calculate P.E. for particle w.r.t. all others
	    for(int j=i+1;j<particles.length;j++){
		
		//calculates separation from i to j

		Vector3d sep= new Vector3d(seperation(particles[i],particles[j]));
		
		//calculate P.E. from i to j
		potential+=(-G*particles[i].getMass()*particles[j].getMass())/sep.mag();
	    }
	}
	
	return potential;
    }
	

    

    /** Create an array of Particles from a Scanner object of the correct format.
     *
     *
     *
     *@param Scanner object to create Particle3d array from.
     *
     *
     *@return Particle3d array object containing all info from Scanner object.
     */

    public static Particle3d[] particleArray (Scanner info){
	Particle3d[] particles = new Particle3d[info.nextInt()];

	//loop to create an array of appropriate size
	for(int i=0;i<particles.length;i++){
	    particles[i]=new Particle3d();
	}
	
	//loop to assign correct infomation to the array
	for(int i=0;i<particles.length;i++){
	    particles[i]=  readParticle(info);

	}

	return particles;
    }


    
    /** Scan variables from input file to create a particle.
     *
     *
     *@param info a Scanner object that provides mass, position, velocity and label.
     *@return Particle3d object created from scanner object.
     */
    
    
    public static Particle3d readParticle(Scanner info){
	Particle3d a= new Particle3d();

	a.setLabel(info.next());
	a.setMass(info.nextDouble());
	a.setPosition(new Vector3d(info.nextDouble(), info.nextDouble(), info.nextDouble()));
	a.setVelocity(new Vector3d(info.nextDouble(), info.nextDouble(), info.nextDouble()));
	
	return a;
    }
    
    /** Method to write file in VMD format
     *
     *
     *@param particles an array of Particle3d objects containing all particles in the system.
     *@param outfile a PrintWriter object pointing towards the file it is desired to write to.
     *@param currentStep an integer specifying the current time-step in the simulation.
     *
     */

    public static void toVMD(Particle3d[] particles, PrintWriter outfile, double currentStep) {
	
	//print number of particles
	outfile.printf("%d\n" ,particles.length);

	//print current time step
	outfile.printf("Point = %10.5f \n" , currentStep);
	
	//print particle infor
	for (int i=0; i < particles.length; i++) {
       
	    outfile.printf("%s\n", particles[i].toString());
	    
	}
	outfile.flush();
    }
    



    /**Determines the relative seperation of the position of two particles
     *
     *
     *@param centre a Particle3d object to determine the seperation from
     *@param orbit a Particle3d object to determine the seperation to
     *
     *@return a Vector3d object representing the relative seperation
     */ 

    // creates position vector pointing from centre to orbit

  public static  Vector3d seperation(Particle3d centre, Particle3d orbit){

	return new Vector3d( Vector3d.subVector(centre.getPosition(), orbit.getPosition()));

    }



    /**Converts a Vector3d object to a normalized vector in the same direction
     *
     *
     *
     *@param a Vector3d object to be normalized
     *
     *@return  Vector3d object that is normalized vector
     */


    public static Vector3d unitHat(Vector3d a){

	double mag=a.mag();

	return( a.scalarDivide(mag));
	    }

    
    /** Method to calculate the Gravitational force between two particles
     *
     *@param particles[] an array of  Particle3d objects to calculate the gravitational force of
     *
     *@return Vector3d[] array containing the Gravitational force of all particles in the system
     */
    
    public static Vector3d[] forceCalc(Particle3d[] particles){

	//Declare and initialize all objects
	Vector3d[] force= new Vector3d[particles.length];
	for(int i=0; i<particles.length; i++){
	    force[i]= new Vector3d();
	}


	double magsep=0.0;
	Vector3d sep= new Vector3d();
	Vector3d newForce= new Vector3d();


	//Loop to cycle through all particles in array
	
	for (int i = 0 ; i < particles.length ; i ++){ 


	    //Loop to calculate force between current particle and all others
	    
	    for(int j=i+1; j<particles.length; j++){
	      
	    sep = Particle3d.seperation(particles[i],particles[j]);

	    magsep = sep.mag();
	    
	    newForce= unitHat(sep).scalarMultiply(particles[i].getMass()*particles[j].getMass()*-1.0*G);
	    
	    newForce= newForce.scalarDivide(magsep*magsep);

	    //sum the current force with previous and save the opposite force in
	    //correct index to eliminate weight of calculations
	force[i]= Vector3d.addVector(newForce, force[i]);
	force[j]= Vector3d.addVector(newForce.scalarMultiply(-1),force[j]);
	   
	    }
	}
	return  (force);



    }


    /**Method to calculate Aphelion and Perihelion of planetary orbits.
     *
     *
     *
     *@param perihelion[] array of doubles to hold the perihelion values for each planet.
     *@param aphelion[] array of doubles to hold the aphelion values for each planet.
     *@param particles[] array of Particle3d objects.
     *
     */

    public static void eccentric(double perihelion[], double aphelion[], Particle3d particles[]){

	//loop to cycle through each particle in system
	for(int j=0; j<particles.length; j++){
	    double separation=seperation(particles[j], particles[0]).mag();
	    
	    //statements to determine whether value should be saved
	    if(separation<perihelion[j]){
		perihelion[j]= separation;
	    }
	    else if(separation>aphelion[j]){
		aphelion[j]= separation;
	    }
	}

}
       



    
    /**Method to calculate the CoM velocity adjustment.
     *
     *@param particles[] an array of particles for the orbitary system.
     *
     */
    public static void adjustedVelocitys(Particle3d[] particles){
	
	//initialise the total mass and center of mass as zero valued
	double massTotal = 0.0;
        Vector3d vCoM = new Vector3d();
	
	//loop to calculate the total mass of the planets
	for(int l=0; l<particles.length; l++){
	     massTotal += particles[l].getMass();
	}
	
	//loop to calculate velocity of the total momentum, calling it vCoM as a tempory place holder
	for(int l=0; l<particles.length; l++){
	    Vector3d momentum = (particles[l].getVelocity()).scalarMultiply(particles[l].getMass());
	    vCoM = Vector3d.addVector(vCoM, momentum);
	    
	}
	
	//assign new value to the velocity of the center of mass
	vCoM = vCoM.scalarDivide(massTotal);
	
	//loop to correct the initial velocity of hte planets 
	for(int l=0; l<particles.length; l++){
	    
	    particles[l].setVelocity(Vector3d.subVector(particles[l].getVelocity(), vCoM));
    
	}

    }
    
    





/**Method to track the orbits of bodies in the system with respect to the Sun (at the origin).
 *
 *
 *
 *
 *@param posInitial array of Vector3d objects containing the current positions of bodies.
 *@param posNext array of Vector3d objects containing the next positions of the bodies.
 *@param orbits array of doubles indicating the number of orbits undergone by the bodies.
 *
 */
    
    public static void orbitTrac(Vector3d[] posInitial, Vector3d[] posNext, double[] orbits){
			     
	//Loop to perform tracking for all particles in system
	for(int i=0; i<posInitial.length; i++){

	    //Calculate dot product of current and next positions for given body
	    double numerator= Vector3d.dot(posInitial[i], posNext[i]);
	    
	    double argument= numerator/(posInitial[i].mag()*posNext[i].mag());
	    
	    //Calculate angle between initial and next position in radians
	    double orbitRad= Math.acos(argument);
	    
	    //convert angle to ratio of a full orbit and increment onto previous saves
	    orbits[i]+= orbitRad/(2*Math.PI);

	    //set new position to current
	    posInitial[i]=posNext[i];
	}
    
    
    
    }


/**Method to track the orbits of the moon around the Earth
 *
 *
 *
 *
 *@param moonInitial a Vector3d object thats the current position of the moon w.r.t. the Earth
 *@param moonNext a Vector3d object thats the next position of the moon w.r.t. the Earth
 *@param moonOrbits a double to track the number of orbits undergone by the moon
 *
 */
    
    public static double moonTrac(Vector3d moonInitial, Vector3d moonNext, double moonOrbits){
			     

	    //Calculate dot product of current and next positions for given body
	    double numerator= Vector3d.dot(moonInitial, moonNext);
	    
	    double argument= numerator/(moonInitial.mag()*moonNext.mag());
	    
	    //Calculate angle between initial and next position in radians
	    double orbitRad= Math.acos(argument);
	    
	    //convert angle to ratio of a full orbit and save
	    moonOrbits= orbitRad/(2*Math.PI);

	    //set new position to current
	    moonInitial=moonNext;


	    return moonOrbits;
    }
}









    

