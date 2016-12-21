package experiment;

import java.lang.reflect.Method;

/**
 * Created by guoy28 on 10/5/16.
 */
public class TestReflection {
    public static void main(String[] args) {
        try {
            Class c = Class.forName("java.sql.Connection");
            Method[] m = c.getDeclaredMethods();
            for (int i = 0; i < 3; i++) {
                System.out.println(m[i]);
            }
        } catch (Throwable e) {
            System.out.println(e.getMessage());
        }
    }
}
