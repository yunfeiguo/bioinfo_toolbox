package yunfeiImplementAlgs4;

import edu.princeton.cs.algs4.StdOut;

public class Date implements Comparable<Date> {
    private static final int[] DAYS = { 0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };

    private final int month;   // month (between 1 and 12)
    private final int day;     // day   (between 1 and DAYS[month]
    private final int year;    // year

    /**
     * Initializes a new date from the month, day, and year.
     * @param month the month (between 1 and 12)
     * @param day the day (between 1 and 28-31, depending on the month)
     * @param year the year
     * @throws IllegalArgumentException if this date is invalid
     */
    public Date(int month, int day, int year) {
        if (!isValid(month, day, year)) throw new IllegalArgumentException();
        this.month = month;
        this.day = day;
        this.year = year;
        System.out.println("Date(int,int,int)");
    }
    
    public Date(String date) {
        System.out.println("Date(String date)");
        String[] f = date.split("/");
        if (f.length != 3) throw new IllegalArgumentException();
        month = Integer.parseInt(f[0]);
        day = Integer.parseInt(f[1]);
        year = Integer.parseInt(f[2]);
        if (!isValid(month, day, year)) throw new IllegalArgumentException();
    }
    public int month() {
        System.out.println("month");
        return month;
    }
    public int day() {
        System.out.println("day");
        return day;
    }
    public int year() {
        System.out.println("year");
        return year;
    }
    private static boolean isValid(int m, int d, int y) {
        System.out.println("isvalid");
        if (m < 1 || m > 12) return false;
        if (d < 1 || d > DAYS[m]) return false;
        if (m == 2 && d == 29 && !isLeapYear(y)) return false;
        return true;
    }
    private static boolean isLeapYear(int y) {
        System.out.println("isleapyear");
        if (y % 400 == 0) return true;
        if (y % 100 == 0) return false;
        return y % 4 == 0;
    }
    public Date next() {
        System.out.println("next");
        if (isValid(month, day + 1, year)) return new Date(month, day + 1, year);
        if (isValid(month + 1, 1, year)) return new Date(month + 1, 1, year);
        return new Date(1, 1, year + 1);
    }
    public boolean isAfter(Date that) {
        System.out.println("isafter");
        return compareTo(that) > 0;
    }
    public boolean isBefore(Date that) {
        System.out.println("isBefore");
        return compareTo(that) < 0;
    }
    @Override
    public int compareTo(Date that) {
        System.out.println("compareto");
        if (this.year < that.year) return -1;
        if (this.year > that.year) return 1;
        if (this.month > that.month) return 1;
        if (this.month < that.month) return -1;
        if (this.day > that.day) return 1;
        if (this.day < that.day) return -1;
        return 0;
    }
    @Override
    public String toString() {
        System.out.println("toString");
        return this.month() + "/" + this.day() + "/" + this.year();
    }
    @Override
    public boolean equals(Object other) {
        System.out.println("equals");
        if (this == other) return true;
        if (other == null) return false;
        if (this.getClass() != other.getClass()) return false;
        Date that = (Date) other; //cannot do:  other = (Date) other becuase 'other' type is Object
        return (this.month() == that.month() && this.year() == that.year()
                && this.day() == that.day());
    }
    @Override
    public int hashCode() {
        int hash = 17;
        hash = 31*hash + this.month;
        hash = 31*hash + this.day;
        hash = 31*hash + this.year;        
        System.out.println("hashCode");
        return hash;
    }
    public static void main(String[] args) {
        Date today = new Date(2, 25, 2004); //test constructor
        StdOut.println(today); //test toString?
        for (int i = 0; i < 10; i++) {
            today = today.next();
            System.out.println(today);
        }
        StdOut.println(today.isAfter(today.next()));
        StdOut.println(today.isAfter(today));
        StdOut.println(today.next().isAfter(today));
        
        Date birthday = new Date("10/16/1971");
        StdOut.println(birthday);
        for (int i = 0; i < 10; i++) {
            birthday = birthday.next();
            StdOut.println(birthday);
        }
    }

}
