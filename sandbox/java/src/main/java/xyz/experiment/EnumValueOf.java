package experiment;

/**
 * Created by guoy28 on 12/15/16.
 */
public class EnumValueOf {
  enum Fruit {
    APPLE, ORGANGE, PEAR;
    @Override
    public String toString() {
      return "Fruit: " + this.name();
    }
  }

  public static void main(String[] args) {
    System.out.println(Fruit.valueOf("APPLE"));
    Food.Apple.serve();
  }

  enum Food {
    Apple {
      public void serve() {
        System.out.println("served after peeling");
      }
    },
    Strawberry {
      public void serve() {
        System.out.println("served after washing");
      }
    };

    public abstract void serve();
    }

}
