package experiment;

import java.util.Arrays;
import java.util.*;
import java.util.OptionalDouble;
import java.util.stream.*;

/**
 * Created by guoy28 on 12/13/16.
 */
public class Java8Map {
  private static class User {
    public int age;
    public String name;
    public String last;
    public int id;
    public User(int age, String f, String l, int i) {
      this.age = age;
      this.name = f;
      this.last = l;
      this.id = i;
    }
  }

  private static List<User> users = Arrays.asList(
          new User(1, "Steve", "Vai", 40),
          new User(4, "Joe", "Smith", 32),
          new User(3, "Steve", "Johnson", 57),
          new User(9, "Mike", "Stevens", 18),
          new User(10, "George", "Armstrong", 24),
          new User(2, "Jim", "Smith", 40),
          new User(8, "Chuck", "Schneider", 34),
          new User(5, "Jorje", "Gonzales", 22),
          new User(6, "Jane", "Michaels", 47),
          new User(7, "Kim", "Berlie", 60)
  );

  public static void main(String[] args) {
    oldJavaWay();
    newJavaWay();
  }
  private static void oldJavaWay() {
    int total = 0;

    for (User u : users) {
      total += u.age;
    }

    double average = (double)total / (double)users.size();

    System.out.println("OLDWAY Average User Age: " + average);
  }

  private static void newJavaWay() {
    double avg = users.parallelStream().mapToInt(u -> u.age).average().getAsDouble();
    System.out.println("newway average: " + avg);
  }
}
