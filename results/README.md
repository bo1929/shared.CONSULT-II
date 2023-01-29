**Format of predictions of CONSULT:**
```
READ_ID:FORM TAX_ID:DIST
```

**Naming convention for CONSULT-predictions:**
- *th*: Threshold for making predictions bottom up, threshold with respect to total vote at root.
- *d1*: Each vote is weighted by $(1-\frac{d}{k})^k$.
