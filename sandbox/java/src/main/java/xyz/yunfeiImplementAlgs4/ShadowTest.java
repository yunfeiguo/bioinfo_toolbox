package yunfeiImplementAlgs4;

public class ShadowTest {
    public int x = 0;
    private class InnerClass {
        private int x = 3;
        private void printX(int x) {
            System.out.println("x: " + x);
            System.out.println("this.x: " + this.x);
            System.out.println("ShadowTest.this.x: " + ShadowTest.this.x);
        }
    }
    public static void main (String[] args) {
        ShadowTest test = new ShadowTest();
        ShadowTest.InnerClass testInner = test.new InnerClass();
        testInner.printX(5);
    }
}

