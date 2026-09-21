import java.util.*;

/**
 * Minimal, dependency-free JSON *parser* (DECISIONS_PHASE4.md section 2.1). Every other JSON need in
 * this codebase is output-only, hand-rolled via the {@code kv}/{@code num}/{@code esc} convention
 * (see {@link JsonExporter}) — but the LLM agent has to actually parse arbitrary, real, nested
 * responses (tool-call argument objects, choices arrays, etc.), which is not safe to do with the
 * flat-body regex extraction (`extractStr`) used elsewhere for simple request bodies. This is a
 * standard recursive-descent parser producing plain {@code Map<String,Object>}/{@code List<Object>}/
 * {@code String}/{@code Double}/{@code Boolean}/{@code null} values — no custom value classes, so
 * callers just use ordinary {@code instanceof}/casts.
 */
public class MiniJson {

    public static Object parse(String s) {
        Parser p = new Parser(s);
        p.skipWs();
        Object v = p.parseValue();
        p.skipWs();
        if (p.pos != p.len) throw new IllegalArgumentException("Trailing content at position " + p.pos);
        return v;
    }

    @SuppressWarnings("unchecked")
    public static Map<String, Object> asObject(Object o) {
        if (o instanceof Map) return (Map<String, Object>) o;
        throw new IllegalArgumentException("Expected JSON object, got: " + (o == null ? "null" : o.getClass()));
    }

    @SuppressWarnings("unchecked")
    public static List<Object> asArray(Object o) {
        if (o instanceof List) return (List<Object>) o;
        throw new IllegalArgumentException("Expected JSON array, got: " + (o == null ? "null" : o.getClass()));
    }

    public static String getStr(Map<String, Object> obj, String key, String dflt) {
        Object v = obj.get(key);
        return v == null ? dflt : String.valueOf(v);
    }

    /**
     * Encodes a plain Map/List/String/Number/Boolean/null value tree back to JSON text — the
     * serialization counterpart to {@link #parse(String)}, used wherever this project's more common
     * hand-rolled-flat-string convention (see {@code JsonExporter}) isn't a good fit because the
     * structure is recursive/variable-shape (e.g. an LLM tool-call's echoed-back arguments, or an
     * agent chat history containing nested tool_calls arrays).
     */
    @SuppressWarnings("unchecked")
    public static String encode(Object v) {
        if (v == null) return "null";
        if (v instanceof String) return encodeString((String) v);
        if (v instanceof Boolean) return String.valueOf(v);
        if (v instanceof Number) {
            double d = ((Number) v).doubleValue();
            if (d == Math.rint(d) && !Double.isInfinite(d) && Math.abs(d) < 1e15) return String.valueOf((long) d);
            return String.valueOf(d);
        }
        if (v instanceof Map) {
            StringBuilder sb = new StringBuilder("{");
            boolean first = true;
            for (Map.Entry<String, Object> e : ((Map<String, Object>) v).entrySet()) {
                if (!first) sb.append(",");
                first = false;
                sb.append(encodeString(e.getKey())).append(":").append(encode(e.getValue()));
            }
            return sb.append("}").toString();
        }
        if (v instanceof List) {
            StringBuilder sb = new StringBuilder("[");
            List<Object> list = (List<Object>) v;
            for (int i = 0; i < list.size(); i++) {
                if (i > 0) sb.append(",");
                sb.append(encode(list.get(i)));
            }
            return sb.append("]").toString();
        }
        throw new IllegalArgumentException("Cannot encode value of type " + v.getClass());
    }

    private static String encodeString(String s) {
        StringBuilder sb = new StringBuilder("\"");
        for (int i = 0; i < s.length(); i++) {
            char c = s.charAt(i);
            switch (c) {
                case '"':  sb.append("\\\""); break;
                case '\\': sb.append("\\\\"); break;
                case '\n': sb.append("\\n"); break;
                case '\r': sb.append("\\r"); break;
                case '\t': sb.append("\\t"); break;
                default:
                    if (c < 0x20) sb.append(String.format("\\u%04x", (int) c));
                    else sb.append(c);
            }
        }
        return sb.append("\"").toString();
    }

    private static class Parser {
        final String s;
        final int len;
        int pos = 0;

        Parser(String s) { this.s = s; this.len = s.length(); }

        void skipWs() {
            while (pos < len && Character.isWhitespace(s.charAt(pos))) pos++;
        }

        char peek() {
            if (pos >= len) throw new IllegalArgumentException("Unexpected end of JSON at position " + pos);
            return s.charAt(pos);
        }

        Object parseValue() {
            skipWs();
            char c = peek();
            switch (c) {
                case '{': return parseObject();
                case '[': return parseArray();
                case '"': return parseString();
                case 't':
                    expectLiteral("true"); return Boolean.TRUE;
                case 'f':
                    expectLiteral("false"); return Boolean.FALSE;
                case 'n':
                    expectLiteral("null"); return null;
                default:
                    return parseNumber();
            }
        }

        void expectLiteral(String lit) {
            if (pos + lit.length() > len || !s.regionMatches(pos, lit, 0, lit.length()))
                throw new IllegalArgumentException("Invalid literal at position " + pos);
            pos += lit.length();
        }

        Map<String, Object> parseObject() {
            Map<String, Object> m = new LinkedHashMap<>();
            pos++; // '{'
            skipWs();
            if (pos < len && s.charAt(pos) == '}') { pos++; return m; }
            while (true) {
                skipWs();
                if (peek() != '"') throw new IllegalArgumentException("Expected string key at position " + pos);
                String key = parseString();
                skipWs();
                if (peek() != ':') throw new IllegalArgumentException("Expected ':' at position " + pos);
                pos++;
                Object val = parseValue();
                m.put(key, val);
                skipWs();
                char c = peek();
                if (c == ',') { pos++; continue; }
                if (c == '}') { pos++; break; }
                throw new IllegalArgumentException("Expected ',' or '}' at position " + pos);
            }
            return m;
        }

        List<Object> parseArray() {
            List<Object> list = new ArrayList<>();
            pos++; // '['
            skipWs();
            if (pos < len && s.charAt(pos) == ']') { pos++; return list; }
            while (true) {
                Object val = parseValue();
                list.add(val);
                skipWs();
                char c = peek();
                if (c == ',') { pos++; continue; }
                if (c == ']') { pos++; break; }
                throw new IllegalArgumentException("Expected ',' or ']' at position " + pos);
            }
            return list;
        }

        String parseString() {
            StringBuilder sb = new StringBuilder();
            pos++; // opening quote
            while (true) {
                if (pos >= len) throw new IllegalArgumentException("Unterminated string");
                char c = s.charAt(pos++);
                if (c == '"') break;
                if (c == '\\') {
                    if (pos >= len) throw new IllegalArgumentException("Unterminated escape");
                    char e = s.charAt(pos++);
                    switch (e) {
                        case '"':  sb.append('"'); break;
                        case '\\': sb.append('\\'); break;
                        case '/':  sb.append('/'); break;
                        case 'b':  sb.append('\b'); break;
                        case 'f':  sb.append('\f'); break;
                        case 'n':  sb.append('\n'); break;
                        case 'r':  sb.append('\r'); break;
                        case 't':  sb.append('\t'); break;
                        case 'u':
                            if (pos + 4 > len) throw new IllegalArgumentException("Bad unicode escape");
                            String hex = s.substring(pos, pos + 4);
                            sb.append((char) Integer.parseInt(hex, 16));
                            pos += 4;
                            break;
                        default:
                            throw new IllegalArgumentException("Bad escape '\\" + e + "' at position " + pos);
                    }
                } else {
                    sb.append(c);
                }
            }
            return sb.toString();
        }

        Object parseNumber() {
            int start = pos;
            if (pos < len && s.charAt(pos) == '-') pos++;
            while (pos < len && Character.isDigit(s.charAt(pos))) pos++;
            if (pos < len && s.charAt(pos) == '.') {
                pos++;
                while (pos < len && Character.isDigit(s.charAt(pos))) pos++;
            }
            if (pos < len && (s.charAt(pos) == 'e' || s.charAt(pos) == 'E')) {
                pos++;
                if (pos < len && (s.charAt(pos) == '+' || s.charAt(pos) == '-')) pos++;
                while (pos < len && Character.isDigit(s.charAt(pos))) pos++;
            }
            if (pos == start) throw new IllegalArgumentException("Invalid number at position " + pos);
            return Double.parseDouble(s.substring(start, pos));
        }
    }
}
