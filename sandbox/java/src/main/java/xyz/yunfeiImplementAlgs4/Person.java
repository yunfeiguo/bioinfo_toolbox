package yunfeiImplementAlgs4;

import edu.princeton.cs.algs4.StdOut;

/******************************************************************************
 *  Compilation:  javac Person.java
 *  Execution:    java Person
 *  Dependencies: StdOut.java
 * 
 *  Implementing equals() in a client-defined type.
 *
 ******************************************************************************/

//public makes this class visible outside the package
//final class makes sure that no one can override the methods defined in this class
//if you use final on a method, then no one can override that specific method.
public final class Person {
    private final String name;
    private final long info;
    
    public Person(String name, long info) {
        this.name = name;
        this.info = info;
    }
    
    //how you're supposed to implement equals
    public boolean equals(Object other) {
        if (other == this) {
            return true;
        }
        if (other == null) {
            return false;
        }
        if (!other.getClass().equals(this.getClass())) {
            return false;
        }
        Person that = (Person) other;
        /*if (!this.name.equals(that.name))
            return false;
        if (this.info != that.info) 
            return false;
        */
        return (this.name.equals(that.name) && this.info == that.info);
    }
    public String toString() {
        return name + " " + info;
    }
    public static void main(String[] args) {
        Person a = new Person("Alice", 1234);
        Person b = new Person("Alice", 1234);
        Person c = new Person("Bob",   1234);
        Person d = new Person("Alice", 4321);
        StdOut.println("a = " + a);
        StdOut.println("b = " + b);
        StdOut.println("c = " + c);
        StdOut.println("d = " + d);
        StdOut.println("a == a: " + a.equals(a));
        StdOut.println("a == b: " + a.equals(b));
        StdOut.println("a == c: " + a.equals(c));
        StdOut.println("a == d: " + a.equals(d));
    }
}
