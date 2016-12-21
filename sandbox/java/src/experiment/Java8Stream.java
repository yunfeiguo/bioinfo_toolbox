package experiment;

import java.util.*;
import java.util.function.Supplier;
import java.util.stream.Collector;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * Created by guoy28 on 12/13/16.
 */
public class Java8Stream {
  static class Person {
    String name;
    int age;

    Person(String name, int age) {
      this.name = name;
      this.age = age;
    }

    @Override
    public String toString() {
      return name;
    }
  }

  static class Foo {
    String name;
    List<Bar> bars = new ArrayList<>();

    Foo(String name) {
      this.name = name;
    }
  }

  static class Bar {
    String name;

    Bar(String name) {
      this.name = name;
    }
  }

  public static void main(String[] args) {
    List<String> l = Arrays.asList("a1", "c1", "bx", "cc");
    l.stream().filter(s -> s.startsWith("c")).
            map(String::toUpperCase).
            sorted().forEach(System.out::println);


    //we don't have to create stream from collections
    Stream.of("a1", "a2", "a3").findFirst().
    ifPresent(System.out::println);

    IntStream.range(0,10).forEach(System.out::print);
    System.out.println();

    Arrays.stream(new int[]{1,2,3,4,5}).
            average().ifPresent(System.out::println);

    Stream.of("a1", "a2","a3").
            map(s -> s.substring(1)).
            mapToInt(Integer::parseInt).
            max().ifPresent(System.out::println);

    IntStream.range(0,10).
            mapToObj(i -> "a" + i).
            forEach(System.out::println);

    Stream.of(1.0, 2.0, 3.0).mapToInt(Double::intValue).
          mapToObj(i -> "double2int" + i).
          forEach(System.out::println);


    Arrays.asList(1,2,3).stream().mapToInt(i->i).average().ifPresent(System.out::println);

    Stream.of(1,2,3).filter(i -> {
      System.out.println("filter: " + i);
      return i < 3;}).forEach(System.out::println);

    Stream.of("a","b","c").map(s -> {
      System.out.println("map: " + s);
      return s.toUpperCase();
    }).filter(s -> {
      System.out.println("filter: " + s);
      return s.startsWith("A");}).forEach(System.out::println);


    Supplier<IntStream> s = () -> Stream.of(1,2,3).mapToInt(i -> i);
    s.get().average().ifPresent(System.out::println);
    System.out.println(s.get().sum());


    List<Person> persons =
            Arrays.asList(
                    new Person("Max", 18),
                    new Person("Peter", 23),
                    new Person("Pamela", 23),
                    new Person("David", 12));

    List<Person> filtered = persons.stream().filter(p -> p.toString().startsWith("M")).collect(Collectors.toList());
    System.out.println(filtered);

    Map<Integer, List<Person>> ageGroup = persons.stream().collect(Collectors.groupingBy(p -> p.age));
    ageGroup.forEach((age, p) -> System.out.format("age %s, %s\n", age, p));

    System.out.println(persons.stream().collect(Collectors.averagingInt(p -> p.age)));


    System.out.println(persons.stream().collect(Collectors.summarizingInt(p -> p.age)));


    System.out.println(persons.stream().map(p -> p.name).collect(Collectors.joining(" and ", "In Germany ", " are here")));

    Map<Integer, String> map = persons.stream().collect(Collectors.toMap(p -> p.age, p -> p.name, (n1, n2) -> n1 + ";" + n2));
    System.out.println(map);

    Collector<Person, StringJoiner, String> personNameCollector =
            Collector.of(
                    () -> new StringJoiner(" | "), //supplier
                    (j, p) -> j.add(p.name.toUpperCase()), //accumulator
                    (j1, j2) -> j1.merge(j2), //combiner
                    StringJoiner::toString); //finisher
    String names = persons.stream().collect(personNameCollector);

    IntStream.range(1, 4)
            .mapToObj(i -> new Foo("Foo" + i))
            .peek(f -> IntStream.range(1, 4)
                    .mapToObj(i -> new Bar("Bar" + i + " <- " + f.name))
                    .forEach(f.bars::add))
            .flatMap(f -> f.bars.stream())
            .forEach(b -> System.out.println(b.name));

    persons.stream().reduce((p1, p2) -> p1.age > p2.age ? p1 : p2).ifPresent(System.out::println);

    Person result = persons.stream().reduce(new Person("", 0), (p1, p2) -> {p1.age += p2.age; p1.name += p2.name; return p1;});
    System.out.println(result);

    System.out.println(persons.stream().reduce(0, (sum, p) -> {sum += p.age;
      System.out.println("accumulator " + sum + " += " + p.age);
      return sum;}, (sum1, sum2) -> {
      System.out.println("combiner " + sum1 + " += " + sum2);
      sum1+=sum2;
      return sum1;
    } ));

    System.out.println(persons.parallelStream().reduce(0, (sum, p) -> {sum += p.age;
      System.out.println("accumulator " + sum + " += " + p.age);
      return sum;}, (sum1, sum2) -> {
      System.out.println("combiner " + sum1 + " += " + sum2);
      sum1+=sum2;
      return sum1;
    } ));
  }
}
