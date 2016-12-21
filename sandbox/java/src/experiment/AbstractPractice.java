package experiment;

/**
 * Created by guoy28 on 12/2/16.
 */
public class AbstractPractice {
  abstract static class MyFirstAbstractClass {
    public MyFirstAbstractClass() {
      System.out.println("abstract class instantiated!?");
    }

    public void say() {
    }
  }

  public static void main(String[] args) {
    MyFirstAbstractClass x = new MyFirstAbstractClass(){
        public void say() {
          System.out.println("say something");
        }
    };
    x.say();
  }
}
