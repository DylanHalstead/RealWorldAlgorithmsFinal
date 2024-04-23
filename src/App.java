import bridges.connect.Bridges;

public class App {
    public static void main(String[] args) throws Exception {
        Bridges bridges = new Bridges(Integer.parseInt(args[0]), args[1], args[2]);
    }
}
