package world.util;

import lombok.AllArgsConstructor;
import lombok.NoArgsConstructor;

/**
 * When using pointers in Java, need to be held in an instance of the object
 * @param <T> Type of Item
 */
@NoArgsConstructor
@AllArgsConstructor
public class Holder<T> {
    private T item;

    public T get() {
        return item;
    }

    public void set(T item) {
        this.item = item;
    }
}
