**Format of predictions of CONSULT:**
```
READ_ID:FORM TAX_ID:DIST
```

**Naming convention for CONSULT-predictions:**
- *th*: Threshold for making predictions bottom up, threshold total vote at root times *th*.
- *d1*: Each vote is weighted by $(1-\frac{d}{k})^k$.
- *c*: Don't make any prediction if the vote threshold is smaller than *c*.
