
import java.io.*;
import java.nio.file.*;
import java.security.*;
import java.util.*;

/**
 * SHA-256 content hashing for cache invalidation.
 * Each base-step manifest records the hash of its inputs;
 * the step is skipped when all input hashes match.
 */
public class ContentHasher {

    public static String hashFile(File file) throws IOException {
        if (!file.exists()) return "MISSING";
        try {
            MessageDigest md = MessageDigest.getInstance("SHA-256");
            byte[] buf = new byte[8192];
            try (InputStream is = new FileInputStream(file)) {
                int n;
                while ((n = is.read(buf)) > 0) md.update(buf, 0, n);
            }
            return hexString(md.digest());
        } catch (NoSuchAlgorithmException e) {
            throw new RuntimeException(e);
        }
    }

    public static String hashString(String value) {
        try {
            MessageDigest md = MessageDigest.getInstance("SHA-256");
            md.update(value.getBytes("UTF-8"));
            return hexString(md.digest());
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    public static String combineHashes(String... hashes) {
        return hashString(String.join("|", hashes));
    }

    public static String combineHashes(List<String> hashes) {
        return hashString(String.join("|", hashes));
    }

    private static String hexString(byte[] bytes) {
        StringBuilder sb = new StringBuilder(bytes.length * 2);
        for (byte b : bytes) sb.append(String.format("%02x", b & 0xff));
        return sb.toString();
    }
}
