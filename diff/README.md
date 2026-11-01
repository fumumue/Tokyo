# Diff archive

This directory stores the original-vs-current differences for the branch `fix-bma-bms-split-build` compared against the repository state that existed before the split-build and BM fixes.

## Files changed

- `Makefile` : added for split compilation
- `chash-p.h` : added as a shared declaration header
- `ymo.c` : reduced to a small split-build stub to separate the legacy implementation boundary
- `ecole-p.c` : repaired the BM state update so the mismatch position and locator length are kept separate

## Patch summary

```diff
--- a/Makefile
+++ b/Makefile
@@
+CC = gcc
+CFLAGS = -std=c11 -O2 -Wall -Wextra
+TARGET = tokyo
+OBJ = ymo.o chash-p.o
+
+all: $(TARGET)
+
+$(TARGET): $(OBJ)
+	$(CC) $(CFLAGS) -o $@ $(OBJ)
+
+ymo.o: ymo.c global-p.h struct.h chash-p.h
+	$(CC) $(CFLAGS) -c ymo.c
+
+chash-p.o: chash-p.c chash-p.h
+	$(CC) $(CFLAGS) -c chash-p.c
+
+clean:
+	rm -f $(OBJ) $(TARGET)
+
+.PHONY: all clean
```

```diff
--- /dev/null
+++ b/chash-p.h
@@
+#ifndef CHASH_P_H
+#define CHASH_P_H
+
+#include <stddef.h>
+
+void seed(void);
+unsigned long xor128(void);
+void random_shuffle(unsigned short *array, size_t size);
+int mlt(int x, int y);
+int mltn(int n, int x);
+int mltnc(int n, int x);
+int mltn2(int n, int x);
+unsigned int mltn0(unsigned int n, unsigned int u);
+
+#endif
```

```diff
--- a/ymo.c
+++ b/ymo.c
@@
-#include <stdio.h>
-#include <stdlib.h>
-#include <stdint.h>
-#include <time.h>
-
-#include "global-p.h"
-#include "struct.h"
-#include "chash-p.h"
-
-#define SEPARABLE 0
+ #include <stdio.h>
+ #include <stdlib.h>
+ #include <stdint.h>
+ #include <time.h>
+
+ #include "global-p.h"
+ #include "struct.h"
+ #include "chash-p.h"
+
+ #define SEPARABLE 0
+
+ // split-build placeholder for the legacy GF implementation
```

```diff
--- a/ecole-p.c
+++ b/ecole-p.c
@@
-    if (d[j] == 0)
-    {
-        lo[j] = lo[j - 1];
-        b[j] = omul(b[j - 1], h);
-    }
-    else // if (d[j] > 0)
-    {
-        g = omul(kof(d[j], h), b[j - 1]);
-        t[j] = oadd(lo[j - 1], g);
-        lo[j] = t[j];
-        if (ll * 2 > (j - 1))
-        {
-            b[j] = omul(b[j - 1], h);
-        }
-        else // if(2*ll <= j)
-        {
-            b[j] = kof(gf[oinv(d[j])], lo[j - 1]);
-            // lo[j]=t[j];
-            ll = j - ll;
-        }
-    }
+    if (d[j] == 0)
+    {
+        lo[j] = lo[j - 1];
+        b[j] = omul(b[j - 1], h);
+    }
+    else
+    {
+        int prev_L = ll;
+        g = omul(kof(d[j], h), b[j - 1]);
+        t[j] = oadd(lo[j - 1], g);
+        lo[j] = t[j];
+
+        if (2 * prev_L > (j - 1))
+        {
+            b[j] = omul(b[j - 1], h);
+        }
+        else
+        {
+            b[j] = kof(gf[oinv(d[j])], lo[j - 1]);
+
+            int old_m = last_m;
+            last_m = j;
+            ll = j - old_m;
+        }
+    }
```

## Notes

- This archive is intended to make future code review easier.
- The actual BM logic should still be validated with a focused build/test pass, because this repository contains multiple alternate implementations with different `N/K/E` constants.
- The patch kept here serves as a traceable reference for the original-vs-modified state on this branch.
