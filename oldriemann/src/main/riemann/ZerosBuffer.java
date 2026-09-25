package riemann;

import java.util.Arrays;

public final class ZerosBuffer {
    private final static double[][] buffer = new double[3][3];
    private static int head = 0;
    private static int tail = 0;
    private static int size = 0;
    private static final int capacity = 3;

    private ZerosBuffer() {
    }

    /**
     * Copies all current elements into the destination array, 
     * ordered from oldest (index 0) to newest.
     * 
     * @param dest The pre-allocated array to copy data into.
     * @return The number of elements actually copied.
     */
    public static int getAllOrdered(double[][] dest) {
        if (dest == null) {
            return 0;
        }
        
        // Copy only up to the available data or destination size to prevent out-of-bounds
        int countToCopy = Math.min(size, dest.length);
        
        for (int i = 0; i < countToCopy; i++) {
            // Calculate the internal index starting from the oldest item (head)
            int internalIndex = (head + i) % capacity;
            dest[i] = buffer[internalIndex];
        }
        
        return countToCopy;
    }

    public static void printAll() {
        
        for (int i = 0; i < size; i++) {
            // Calculate the internal index starting from the oldest item (head)
            int internalIndex = (head + i) % capacity;
            System.out.println(i + " " +Arrays.toString(buffer[internalIndex]));
            //buffer[internalIndex];
        }
        
    }
    
    public static void reset() {
    	head = tail = size = 0;
    }

    public static double[] getRow(int index) {
    	if (index > size -1) {
    		throw new IllegalArgumentException();
    	}
        int internalIndex = (head + index) % capacity;
        return buffer[internalIndex];
    	
    }
    
    public static void put(double[] value) {
        buffer[tail] = value;
        tail = (tail + 1) % capacity;

        if (size == capacity) {
            head = (head + 1) % capacity; 
        } else {
            size++;
        }
    }

    public static boolean isEmpty() {
        return size == 0;
    }

    public static int size() {
        return size;
    }
}
