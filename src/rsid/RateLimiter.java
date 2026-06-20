package rsid;

/**
 * Simple token-bucket rate limiter. Thread-safe.
 */
public class RateLimiter {
    private final double tokensPerSec;
    private double tokens;
    private long lastRefill;

    public RateLimiter(double requestsPerSecond) {
        this.tokensPerSec = requestsPerSecond;
        this.tokens = requestsPerSecond;
        this.lastRefill = System.nanoTime();
    }

    public synchronized void acquire() throws InterruptedException {
        refill();
        while (tokens < 1.0) {
            long waitMs = (long) ((1.0 - tokens) / tokensPerSec * 1000) + 1;
            Thread.sleep(waitMs);
            refill();
        }
        tokens -= 1.0;
    }

    private void refill() {
        long now = System.nanoTime();
        double elapsed = (now - lastRefill) / 1_000_000_000.0;
        tokens = Math.min(tokensPerSec, tokens + elapsed * tokensPerSec);
        lastRefill = now;
    }
}
