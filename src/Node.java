public class Node {
    private double x;
    private double y;
    private int order;
    public double getX() {
        return x;
    }

    public double getY() {
        return y;
    }

    public int getOrder() {
        return order;
    }

    public void setX(double x) {
        this.x = x;
    }

    public void setY(double y) {
        this.y = y;
    }

    public void setOrder(int order) {
        this.order = order;
    }

    Node(double x, double y,int order){
        this.x=x;
        this.y=y;
        this.order=order;
    }
    public double distanceTo(Node other) {
        double dx = x - other.x;
        double dy = y - other.y;
        return Math.sqrt(dx * dx + dy * dy);
    }

}
