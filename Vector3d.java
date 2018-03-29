/** 
 *A class for 3-D Vectors, complete with constructors, setters and getters,
 *methods for two vector, subtraction, multiplication and
 * division by a scalar, scalar product, vector product, magnitude, 
 *and magnitude squared as well as a toString() method
 *
 *@author C.W.John
 *@author M.O'Shea
 */



public class Vector3d {


    private double xx;
    private double yy;
    private double zz;


    /**Constructor for  empty vector if object called with no arguments.
    */
    public Vector3d() {
	this.setVector(0.0, 0.0, 0.0);
    }

    /**Constructor for copying the vector presented as the argument to a new vector.
     *@param orig vector to copy
    */
    public Vector3d(Vector3d orig){
	this.setVector(orig.getX(), orig.getY(), orig.getZ());
    }

    /**Constructor for explicitly defining components of vector with the arguments.
     *@param x to set value of vector
     *@param y to set value of vector
     *@param z to set value of vector
    */
    public Vector3d(double x,double y,double z){
	this.setVector(x,y,z);
    }


    /**Setter to assign internal x-variable of a vector object.
     *@param x to set value of component
    */

    public void setX(double x){
	this.xx= x;
    }

    /**Setter to assign internal y-variable of a vector object.
     *@param y to set value of component
    */
    public void setY(double y){
	this.yy=y;
    }

    /**Setter to assign internal z-variable of a vector object.
     *@param z to set value of component
    */
    public void setZ(double z){
	this.zz=z;
    }

    /**Getter to retrieve internal x-variable of vector object.
     *@return a double of internal x value
     */

    public double getX(){
	return this.xx;
    }

     /**Getter to retrieve internal y-variable of vector object.
      * @return a double of internal y value
     */
    public double getY(){
	return this.yy;
    }

     /**Getter to retrieve internal z-variable of vector object.
      * @return a double of internal z value
     */
    public double getZ(){
	return this.zz;
    }

    /** Setter to create an entire 3d-vector.
     *@param x to set value of component
     *@param y to set value of component
     *@param z to set value of component
    */

    public void setVector(double x, double y, double z){
	this.setX(x);
	this.setY(y);
	this.setZ(z);
    }

    
    /* instance methods for magnitude squared, magnitude of vector, and prints a string */
    
    
    /**Instance method to return magnitude squared of a vector.
     *@return double value of magnitude squared
     */
    public double magSquared(){
	return this.getX()*this.getX() + this.getY()*this.getY()+this.getZ()*this.getZ();
    }
    
    /**Instance method to return magnitude of a vector.
     *@return double value of magnitude squared
     */
    public double mag(){
	return Math.sqrt(this.magSquared());
    }


    /**Instance method to convert vector into a readable string.
     *Methods called 'toString' are automatically called when
     *an object is printed. <br>
     *@return String vector formatted to string
     *
     */
      public String toString(){
	double x= this.getX();
	double y= this.getY();
	double z= this.getZ();


	return Double.toString(x)+ " " + Double.toString(y) +" " + Double.toString(z);
	    }



/* Instance methods for scalar multiplication and scalar division by a double
    */

    /**Instance method for multiplication of a vector with a scalar value.
     *@param a double to multiply vector by
     */
    public Vector3d scalarMultiply(double a){
       

    return new Vector3d(this.getX()*a, this.getY()*a, this.getZ()*a);
	
    }


    /**Instance method for divsion of a vector with a scalar value.
     *@param a double to divide vector by
     */

    public Vector3d scalarDivide(double a){
       

    return new Vector3d(this.getX()/a, this.getY()/a, this.getZ()/a);
	
    }
    
    /**Instance method for negating a vector.
     */

    public Vector3d negate(){
	
	return new Vector3d(xx*-1.0, yy*-1.0, zz*-1.0);
    }

    public Vector3d unit() {
	return scalarMultiply(1.0/mag());
    }
    



    /* Static methods to perform vector addition, subtraction, cross product, dot product
	on two 3d-Vectors */

    public void add(Vector3d v) {
	xx += v.getX();
	yy += v.getY();
	zz += v.getZ();
    }

    /** Static method to perform vector addition of two 3D vectors.
     *@param Vector3d object to add
     *@param Vector3d object to add
     *@return Vector3d object of sum of vectors
 */

    public static Vector3d addVector(Vector3d a, Vector3d b){
	return new Vector3d(a.getX()+ b.getX(), a.getY()+b.getY(), a.getZ()+b.getZ());
    }

    /** Static method to perform vector subtraction of two 3D vectors.
     *@param Vector3d object to substract
     *@param Vector3d object to subtract
     *@return Vector3d object of difference between vectors
 */

    public static Vector3d subVector(Vector3d a, Vector3d b){
	return new Vector3d(a.getX()- b.getX(), a.getY()-b.getY(), a.getZ()-b.getZ());
    }
	

    /** Static method to perform cross product of two 3D vectors.
     *@param Vector3d object to cross
     *@param Vector3d object to cross
     *@return Vector object of cross product of parameters
     */

    public static Vector3d cross(Vector3d a, Vector3d b){
	return new Vector3d((a.getY()*b.getZ() - b.getY()*a.getZ()), (b.getX()*a.getZ()- a.getX()*b.getZ()), (a.getX()*b.getY()- b.getX()*a.getY()));
    }


     /** Static method to perform cross product of two 3D vectors. 
      *@param Vector object to dot product
      *@param Vector object to dot product
      *@return Double value of dot product between vectors
      */

    public static double dot(Vector3d a, Vector3d b){
	return ((a.getX()*b.getX())+(a.getY()*b.getY())+(a.getZ()*b.getZ()));
    }

};

