# Maaslin2_volcano_plots

### Table of Contents

<aside>

# Part 1. Volcano Plots Using Discriminant Taxa

</aside>

### 1. Loading Libraries

```r
library(ggplot2)
library(ggrepel)
library(dplyr)
library(readxl)
```

These lines load the necessary libraries:

- `ggplot2` for creating plots.
- `ggrepel` for better label placement in plots.
- `dplyr` for data manipulation.
- `readxl` for reading Excel files.

### **2. Reading Data**

```r
# Read in main dataset
data <- read_excel("Maaslin_allresults_COMBINED.xlsx", sheet = "SF v TP") %>% as.data.frame()

# Read in discriminant taxa from Maaslin_discrim_results.xlsx
discriminant_otus <- read_excel("Maaslin_discrim_results.xlsx", sheet = "SF5_discrim") %>%
  select(OTUID) %>% 
  pull()
```

- The first block reads the main dataset from the specified Excel file and sheet, converting it to a data frame.
- The second block reads a list of discriminant taxa from another Excel file and sheet, selecting only the `OTUID` column and converting it to a vector.
  - **Why Convert to Data Frame and Vector?**

    | **Conversion**              | **Reason**                                                     |
    | --------------------------------- | -------------------------------------------------------------------- |
    | **Tibble → Data Frame**    | Ensures compatibility with standard R functions                      |
    | **Single Column → Vector** | Allows easy filtering, comparisons, and membership checks `(%in%)` |

    When reading data from Excel using read_excel(), the resulting object is typically a **tibble**, which is a type of data frame with some additional features. In many cases, we convert it to a **regular data frame** or **vector** for easier manipulation.


    1. ***Converting the Main Dataset to a Data Frame***

       `data <- read_excel("Maaslin_allresults_COMBINED.xlsx", sheet = "SF v TP") %>% as.data.frame()`

       **A. Ensures compatibility** with functions that expect a standard data.frame instead of a tibble.

       **B.** Some older R functions (e.g., dplyr::filter(), ggplot2) may have issues with tibbles, so converting ensures no unexpected behavior.

       **C. Standardizes the structure** if merging or manipulating the dataset later.
    2. **Converting the OTUID Column to a Vector**

       `discriminant_otus <- read_excel("Maaslin_discrim_results.xlsx", sheet = "SF5_discrim") %>% select(OTUID) %>% pull()`

       **Why?**

       **A. Extracts only one column (OTUID)** from the dataset.

       **B. Turns it into a vector** instead of keeping it as a data frame (or tibble) with one column.

       C. This allows easier comparisons, such as:

       `data$OTUID %in% discriminant_otus`  → Checks if OTUs in `data` are also in `discriminant_otus`

       D. **More efficient filtering:** When we filter or subset data, it’s easier and faster to work with a vector than a single-column data frame.

### **3. Data Preparation**

```r
# Ensure matching format
data$OTUID <- trimws(as.character(data$OTUID))
discriminant_otus <- trimws(as.character(discriminant_otus))
```

These lines ensure that the `OTUID` values in both datasets are trimmed of whitespace and converted to character strings for consistent comparison.

### **4. Data Transformation**

```r
data <- data %>%
  mutate(log_q = -log10(qval),  # this line can be done in excel or using R
         Enrichment = case_when(
           qval >= 0.05 ~ "Non-Significant",  # Light lavender for non-significant
           coef > 0 ~ "SF Enriched",  # Soft blue for SF Enriched
           coef < 0 ~ "TP Enriched"   # Peach for TP Enriched
         ),
         Highlight = ifelse(qval < 0.05 & OTUID %in% discriminant_otus, "Highlighted", "Regular"))  # Keep only significant consistent OTUs highlighted
```

- `mutate(log_q = -log10(qval))`: Adds a new column `log_q` which is the negative log10 of the `qval` column.

  - **Why Use `mutate(log_q = -log10(qval))`?**

    - **This makes it easier to set a threshold** (e.g., log_q > 1.3 corresponds to qval < 0.05).
    - **The higher the log_q value, the more significant the result.**
    - Used in **volcano plots**, where log_q is plotted on the y-axis to highlight significant points.

    When analyzing statistical significance in biological and microbiome studies, it's common to **transform p-values (or q-values)** using the negative log10 scale. **This transformation makes differences in significance more visually interpretable.**

    ---

    ### **Mathematical Explanation**


    - The **q-value (`qval`)** represents the adjusted p-value (false discovery rate adjusted).
    - Since q-values are typically **very small numbers**, plotting them directly can be **difficult to interpret**.
    - **Taking `log10(qval)` transforms these small numbers into larger, more intuitive values**.

    ### **Example Transformation**

    | **qval (Raw Value)** | **log_q = -log10(qval)** |
    | -------------------------- | ------------------------------ |
    | `0.1`                    | `1.0`                        |
    | `0.01`                   | `2.0`                        |
    | `0.001`                  | `3.0`                        |
    | `0.0001`                 | `4.0`                        |

    - **Larger `log_q` values indicate greater statistical significance**.
    - This transformation makes **highly significant points stand out** in plots.

    ### **Application in Volcano Plots**

    ```r
    data <- data %>%
      mutate(log_q = -log10(qval))  # Transform q-values for better visualization
    ```
- `Enrichment = case_when(...)`: Adds a new column `Enrichment` based on the `qval` and `coef` values, categorizing them as "Non-Significant", "SF Enriched", or "TP Enriched".
- `Highlight = ifelse(...)`: Adds a new column `Highlight` to mark significant and consistent OTUs as "Highlighted".

### **5. Plotting**

```r
# Base volcano plot
plot <- ggplot(data, aes(x = coef, y = log_q, color = Enrichment)) +
  geom_point(size = 1.5, alpha = 0.9) +  # Small but visible dots

  # Cute & colorblind-friendly color mapping
  scale_color_manual(values = c(
    "SF Enriched" = "#CC79A7",  
    "TP Enriched" = "#56B4E9",  
    "Non-Significant" = "gray"  
  )) +

  # Overlay only significant, consistent OTUs with sunny yellow
  geom_point(data = data %>% filter(Highlight == "Highlighted"), 
             aes(x = coef, y = log_q), color = "black", size = 1.5, alpha = .2, shape = 16) +  

  # Dashed vertical line at x = 0, dotted horizontal at significance threshold
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "gray") +

  # Labels & theme adjustments
  labs(x = "Coefficient", y = "-log10(q-value)", color = "Enrichment") +
  theme_minimal() +
  theme(legend.position = "top")

# Show plot
print(plot)
```

- `ggplot(data, aes(x = coef, y = log_q, color = Enrichment))`: Initializes the plot with `coef` on the x-axis, `log_q` on the y-axis, and colors based on `Enrichment`.
- `geom_point(size = 1.5, alpha = 0.9)`: Adds points to the plot with specified size and transparency.
- `scale_color_manual(values = c(...))`: Manually sets the colors for the different enrichment categories.
- `geom_point(data = data %>% filter(Highlight == "Highlighted"), ...)`: Overlays additional points for highlighted OTUs.
- `geom_vline(xintercept = 0, linetype = "dashed", color = "black")`: Adds a dashed vertical line at x = 0.
- `geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "gray")`: Adds a dotted horizontal line at the significance threshold.
- `labs(x = "Coefficient", y = "-log10(q-value)", color = "Enrichment")`: Sets the axis labels and legend title.
- `theme_minimal() + theme(legend.position = "top")`: Applies a minimal theme and positions the legend at the top.
- `print(plot)`: Displays the plot.

This code creates a volcano plot to visualize the enrichment of OTUs, highlighting significant and consistent OTUs.

<aside>

# Part 2. Incorporating Core Taxa

Volcano Plot with Core and Discriminant Taxa

</aside>

## Why Look at a Volcano Plot with Both Core and Discriminant Taxa?

Adding both core taxa and discriminant taxa to the volcano plot provides more ecological and statistical insights into the microbial community.

### 1. Identifies Significant Taxa Driving Differences Between SF and TP

- The volcano plot already shows which taxa are significantly enriched in SF or TP.
- By focusing on q-values (significance) and coefficients (effect size), we see which taxa are statistically important in differentiating the two sites.

**Benefit:** This tells us which taxa are the main drivers of differences between the sites.

### 2. Highlights Taxa That Are Both Discriminant and Core

- Discriminant taxa (black solid dots) are statistically significant in separating SF and TP (based on independent analysis).
- Core taxa (circled points) are found across all sites and depths, meaning they are ecologically consistent.

**Benefit:** If a taxon is both core and discriminant, it is likely a keystone taxon—playing an important role across all sites while also showing enrichment at one site.

### 3. Shows Whether Core Taxa Are Also Driving Site Differences

- Some core taxa may be significantly enriched in SF or TP, while others may be evenly distributed.
- Seeing a core taxon highly enriched in one site suggests it thrives better in that environment.

**Benefit:** Helps understand if core taxa are site-specific indicators or environmentally flexible.

### 4. Helps Interpret Functional Roles of Core Taxa

- Core taxa are often foundational members of the microbiome.
- By labeling them with taxonomy strings, we can infer potential functional roles based on known taxonomy.

**Benefit:** Provides biological meaning to patterns in microbiome structure.

## Insights You Can Get from This Plot

### 1. Which core taxa are significantly enriched in SF or TP?

- If a core taxon is enriched in SF or TP, it may indicate habitat preference despite being widespread.

### 2. Which discriminant taxa are also core?

- If a taxon is both discriminant and core, it may be an important ecological driver.

### 3. Are core taxa mostly non-significant?

- If most core taxa are gray (non-significant), they may not be major drivers of site differences, but rather stable background members of the microbiome.

### 4. Functional Implications of Core Taxa

- Seeing the taxonomy labels helps infer their possible ecological roles.
- For example, if **Nitrospirae** is core and enriched in SF, it might indicate **nitrification differences** between sites.

---

<aside>

# **Updating the Volcano Plot to Use d70 Core Taxa and Add a Legend**

</aside>

*Now, let’s modify the code to label core taxa with their taxonomy strings on the volcano plot → Family.*

### **1. Read and Merge d70 Core Taxa with Family Information**

```r
# Read d70 core taxa and select OTUID and Family column
core_taxa_d70 <- read_excel("CoreTaxa_70_16S_ASCC.xlsx", sheet = "d70") %>%
  select(OTUID, Family)

# Merge taxonomy information into the main dataset
data <- data %>%
  left_join(core_taxa_d70, by = "OTUID") %>%  # Add Family column from d70 core taxa
  mutate(Family_Label = case_when(
    OTUID %in% core_taxa_d70$OTUID ~ Family,  # Assign Family names only to d70 core taxa
    TRUE ~ NA_character_  # Keep all other points unlabeled
  ))

```

---

### **2. Modify the Volcano Plot to Label Core Taxa and Add a Legend**

```r
# Base volcano plot
plot <- ggplot(data, aes(x = coef, y = log_q, color = Enrichment)) +
  geom_point(size = 1.5, alpha = 0.9) +

  # Custom color scheme for enrichment
  scale_color_manual(values = c(
    "SF Enriched" = "#CC79A7",
    "TP Enriched" = "#56B4E9",
    "Non-Significant" = "gray"
  )) +

  # Overlay discriminant taxa (solid red points)
  geom_point(data = data %>% filter(Highlight == "Highlighted"),
             aes(x = coef, y = log_q, shape = "Discriminant Taxa"), size = 3, color = "#AB3329", alpha = 0.8) +

  # Overlay core taxa (circled points)
  geom_point(data = data %>% filter(OTUID %in% core_taxa_d70$OTUID),
             aes(x = coef, y = log_q, shape = "Core Taxa"), fill = NA, size = 2, stroke = 1, color = "black") +

  # Add Family labels for core taxa
  geom_text_repel(data = data %>% filter(!is.na(Family_Label)),
                  aes(x = coef, y = log_q, label = Family_Label),
                  size = 3, color = "black", fontface = "italic", box.padding = 0.5, segment.color = "gray") +

  # Define shape legend
  scale_shape_manual(values = c("Discriminant Taxa" = 16, "Core Taxa" = 21),
                     name = "Taxa Type") +

  # Threshold lines
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "gray") +

  # Labels & theme
  labs(x = "Coefficient", y = "-log10(q-value)", color = "Enrichment") +
  theme_minimal() +
  theme(legend.position = "top")

# Show plot
print(plot)

```

---

## **Explanation of Changes**

- **Now uses d70 core taxa (core across all sites and depths).**
- **Labels core taxa with their Family names** from the d70 core dataset.
- **Distinguishes core and discriminant taxa in the legend:**
  - **Discriminant Taxa** = **Solid points (`shape = 16`)**
  - **Core Taxa** = **Circled points (`shape = 21`)**
- **Uses `geom_text_repel()` to prevent label overlap and ensure clarity.**
- **Adds a legend to clearly show what the different markers mean.**

---

**What’s New in This Version?**

- **Core taxa (d70) are labeled** with their taxonomy string (`asv_string_silva`).
- **Uses `geom_text_repel()`** so labels don’t overlap.
- **Circles core taxa** while keeping solid black dots for discriminant taxa.
- **Maintains color-coded enrichment categories** (SF vs. TP).

## **Final Thoughts**

This update:

- **Ensures core taxa (from d70) are labeled with their Family names.**
- **Keeps the volcano plot clean while providing key taxonomic insights.**
- **Clearly differentiates between core taxa (circled) and discriminant taxa (solid) in the legend.**
- **Prevents text overlap for better readability.**

<aside>

# How Would Using Other Core Sheets Change the Analysis?

</aside>

### 1. d70 (Current Approach) - Core Across All Sites and Depths

- This version highlights taxa that are consistently found across all sites (SF, TP, SJ) and depths.
- **Advantage:** Shows truly universal core taxa that are found everywhere.
- **Limitation:** Some site-specific important taxa might not be labeled.

### 2. d70_SF & d70_TP - Core Taxa Within SF and TP Separately

- If you use **d70_SF**, you only highlight core taxa for **SF**.
- If you use **d70_TP**, you only highlight core taxa for **TP**.
- **Advantage:** Helps show which core taxa are site-specific.
- **Limitation:** Some important taxa might be missed if they are **core in both sites but appear only in one list**.

### 3. d70_SF + d70_TP - Site-Specific Core Taxa

- If you **combine** core taxa from **d70_SF and d70_TP**, you highlight taxa that are **core within at least one of these sites**.
- **Advantage:** Captures core taxa that are stable within SF or TP but not necessarily across all sites.
- **Limitation:** You lose the ability to distinguish whether a taxon is **core to both sites or just one**.

### 4. Depth-Specific Core Taxa (d70_shallow, d70_deep, d70_OM)

- These sheets highlight taxa that are **core within specific depth layers**.
- **Advantage:** Helps separate taxa based on **where they are found in the soil profile**.
- **Limitation:** Might not align well with a **SF vs. TP comparison**, since it’s based on depth rather than site.

## How Would This Affect the Volcano Plot?

- **d70 (universal core taxa)** → Core taxa appear across **all sites and depths** and are **circled** in the plot.
- **d70_SF / d70_TP (site-specific core taxa)** → You could:

  - **Circle SF core taxa differently from TP core taxa.**
  - **Use different shapes for SF vs. TP core taxa.**
- **Depth-specific core taxa** → Could add **depth-based color coding** (e.g., different fill colors for different depths).

  ---

## Modified Code for Different Core Taxa Comparisons

### **1. d70_SF & d70_TP - Site-Specific Core Taxa**

```r
core_taxa_sf <- read_excel("CoreTaxa_70_16S_ASCC.xlsx", sheet = "d70_SF") %>%
  select(OTUID, asv_string_silva)

core_taxa_tp <- read_excel("CoreTaxa_70_16S_ASCC.xlsx", sheet = "d70_TP") %>%
  select(OTUID, asv_string_silva)

```

**Explanation:**

- Reads **core taxa for SF and TP separately**.
- Can be used to highlight **site-specific** taxa.

---

### **2. d70_SF + d70_TP - Combined Site-Specific Core Taxa**

```r
core_taxa_combined <- bind_rows(core_taxa_sf, core_taxa_tp) %>%
  distinct(OTUID, .keep_all = TRUE)

```

**Explanation:**

- Combines **both SF and TP core taxa** into one dataset.
- Ensures **each OTUID appears only once**.

---

### **3. Depth-Specific Core Taxa (d70_shallow, d70_deep, d70_OM)**

```r
core_taxa_shallow <- read_excel("CoreTaxa_70_16S_ASCC.xlsx", sheet = "d70_shallow") %>%
  select(OTUID)

core_taxa_deep <- read_excel("CoreTaxa_70_16S_ASCC.xlsx", sheet = "d70_deep") %>%
  select(OTUID)

core_taxa_om <- read_excel("CoreTaxa_70_16S_ASCC.xlsx", sheet = "d70_OM") %>%
  select(OTUID)

```

**Explanation:**

- Reads **core taxa by depth layer**.
- Could be used to **color-code depths in the volcano plot**.

---

## Changes in This Version

- **Site-specific core taxa (SF & TP) are labeled separately.**
- **Uses different markers:**
  - **SF Core:** Pink circles
  - **TP Core:** Blue squares
- **Adds taxonomy labels for core taxa.**

---

## Insights From This Version

### 1. Are SF and TP core taxa different?

- If SF and TP core taxa separate clearly, this suggests **environment-specific microbiomes**.

### 2. Are core taxa enriched in specific sites?

- If **SF Core taxa** are also SF-enriched, they might play a **key role in SF soils**.

### 3. Which site has more core taxa driving differences?

- If **SF Core taxa** are highly enriched in SF but **TP Core taxa are mostly non-significant**, then SF might have a **more unique microbial community**.

---

## Final Thoughts

This version of the volcano plot lets you:

- **Visualize core taxa specific to each site (SF & TP).**
- **Compare site-specific core taxa with discriminant taxa.**
- **Label and separate taxa by taxonomy and ecological role.**

<aside>

## **Updating the Volcano Plot to Use d70 Core Taxa and Add a Legend (very detailed explanation)**

**Step 1: Load Required Libraries**

```r
library(readxl)  # For reading Excel files
library(dplyr)   # For data manipulation
library(ggplot2) # For plotting
```

- **Why?**
  - **readxl** allows us to import Excel files into R.
  - **dplyr** makes it easy to manipulate and filter data.
  - **ggplot2** is used to create the volcano plot.

**Step 2: Read in the Main Dataset**

```r
data <- read_excel("Maaslin_allresults_COMBINED.xlsx", sheet = "SF v TP") %>% as.data.frame()
```

- **Why?**
  - This loads the **Maaslin2 results** from the “SF v TP” sheet into a **data frame**.
  - This dataset contains key statistical results, including:
  - **OTUID** (unique identifier for taxa)
  - **qval** (p-value adjusted for multiple comparisons)
  - **coef** (effect size: tells us which site a taxon is enriched in)

**Step 3: Read in Discriminant Taxa**

```r
discriminant_otus <- read_excel("Maaslin_discrim_results.xlsx", sheet = "SF5_discrim") %>%
  select(OTUID) %>%
  pull()
```

- **Why?**
  - We **import** a list of **discriminant taxa** from a different analysis (Maaslin2 discriminant results).
  - **select(OTUID) %>% pull()** extracts only the **OTUID column** as a simple character vector.
  - These taxa will be **highlighted in the plot** to show that they were also significant in a separate analysis.

**Step 4: Read in Core Taxa (Universal Across All Sites and Depths)**

```r
core_taxa <- read_excel("CoreTaxa_70_16S_ASCC.xlsx", sheet = "d70") %>%
  select(OTUID) %>%
  pull()
```

- **Why?**
  - The **d70** sheet contains **core taxa across all sites and depths**.
  - We extract **OTUIDs** of these core taxa to check which ones appear in our main dataset.
  - These taxa will be **circled in the plot** to indicate that they are universal.

**Step 5: Ensure Matching Format for OTUIDs**

```r
data$OTUID <- trimws(as.character(data$OTUID))
discriminant_otus <- trimws(as.character(discriminant_otus))
core_taxa <- trimws(as.character(core_taxa))
```

- **Why?**
  - **Ensures all OTUIDs are in the same format** (character strings with no extra spaces).
  - This is crucial because different formatting (e.g., extra spaces, factors vs. characters) can **cause mismatches** when filtering or merging datasets.

**Step 6: Add Metadata to the Data for Plotting**

```r
data <- data %>%
  mutate(log_q = -log10(qval),  # Transform q-value
         Enrichment = case_when(
           qval >= 0.05 ~ "Non-Significant",  # Non-significant taxa (gray)
           coef > 0 ~ "SF Enriched",          # Taxa enriched in SF (pink)
           coef < 0 ~ "TP Enriched"           # Taxa enriched in TP (blue)
         ),
         Highlight = ifelse(qval < 0.05 & OTUID %in% discriminant_otus, "Highlighted", "Regular"),
         Core_Taxa = ifelse(OTUID %in% core_taxa, "Core", "Non-Core"))  # Mark core taxa
```

- **Why?**
  1. **log_q = -log10(qval)**
     1. Converts **q-values** into **-log10 scale** to make small p-values (high significance) more visible.
     2. Higher log_q values = **stronger significance**.
  2. **Assigns “Enrichment” categories**
     1. **qval >= 0.05** → **“Non-Significant”** (gray)
     2. **coef > 0** → **“SF Enriched”** (pink)
     3. **coef < 0** → **“TP Enriched”** (blue)
     4. This **color-codes points** in the plot.
  3. **Highlights Discriminant Taxa**
  4. If an OTU appears in **both** the main dataset and the **discriminant list**, it is **highlighted**.
  5. **Marks Core Taxa**
     1. If an OTU appears in the **d70 core taxa list**, it is **labeled as “Core.”s**
     2. These will be **circled** in the plot.

 **Step 7: Create the Volcano Plot**

```r
plot <- ggplot(data, aes(x = coef, y = log_q, color = Enrichment)) +
  geom_point(size = 1.5, alpha = 0.9) +
```

- **Why?**
  - The **x-axis (coef)** represents effect size (which site the taxon is enriched in).
  - The **y-axis (-log10(q-value))** represents statistical significance.
  - Each point represents a taxon.

**Step 8: Custom Color Scheme**

```r
  scale_color_manual(values = c(
    "SF Enriched" = "#CC79A7",
    "TP Enriched" = "#56B4E9",
    "Non-Significant" = "gray"
  )) +
```

- **Why?**
  - **Custom colors** for **SF (pink), TP (blue), and Non-Significant (gray)** to make the plot more readable.

**Step 9: Overlay Discriminant Taxa (Black Dots)**

```r
  geom_point(data = data %>% filter(Highlight == "Highlighted"), 
             aes(x = coef, y = log_q), color = "black", size = 1.5, alpha = .2, shape = 16) +
```

- **Discriminant taxa** are plotted in **solid black** to indicate they were also significant in another analysis.

**Step 10: Overlay Core Taxa (Circled Points)**

```r
geom_point(data = data %>% filter(Core_Taxa == "Core"),
aes(x = coef, y = log_q), shape = 21, fill = NA, size = 2, stroke = 1, color = "black") +
```

- **Why?**
  - **Core taxa (from d70)** are **circled (open markers, black outline)** to highlight their universal presence across all sites and depths.

**Step 11: Add Threshold Lines**

```r
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "gray") +
```

- **Why?**
  - **Vertical dashed line at x = 0** → Separates SF-enriched and TP-enriched taxa.
  - **Horizontal dotted line at y = -log10(0.05)** → Represents significance threshold (q < 0.05).

**Step 12: Labels and Theme**

```r
  labs(x = "Coefficient", y = "-log10(q-value)", color = "Enrichment") +
  theme_minimal() +
  theme(legend.position = "top")
```

- **Why?**
  - Adds **axis labels** and **a minimal theme** for clarity.
  - Moves **legend to the top** to improve layout.

 **Step 13: Print the Plot**

```r
print(plot)
```

- **Why?**

  - Displays the final volcano plot!
- **Summary of New Features:**

  - **Core taxa** are **circled** (black outline, transparent inside).
  - **Significant discriminant taxa** remain **highlighted** (black solid points).
  - The volcano plot still distinguishes **SF Enriched**, **TP Enriched**, and **Non-Significant** taxa.

  The **core taxa** are identified using the **OTUID** column, which matches your main dataset. To incorporate them into your volcano plot, I’ll:

  - Extract **OTUIDs** from both **SF** and **TP** core taxa sheets.
  - Mark them in your dataset as **“Core Taxa”**.
  - Modify the plot to **circle** these core taxa (using open circles or another visual distinction).

  **Changes:**

  - Reads core taxa from the **d70_SF** and **d70_TP** sheets.
  - Adds a new column "Core_Taxa" to indicate core taxa presence.
  - Uses geom_point(shape = 21, fill = NA, size = 2, stroke = 1) to **circle** core taxa.

<aside>

## Final Thoughts

This volcano plot now provides both statistical significance and ecological meaning:

- **Taxa driving SF vs. TP differences** are **color-coded**.
- **Key core taxa are labeled**, helping understand their roles.
- **Discriminant taxa are highlighted**, connecting statistical and ecological insights.

Let me know if you need any refinements.

</aside>

</aside>

<aside>

## Updated R Code (with taxonomic string labels, optional)

## Adding Taxonomy Labels to Core Taxa (optional)

We will now update the volcano plot to:

- Use **d70 core taxa** (universal across all sites and depths).
- **Label core taxa** using their **Family** taxonomy (optional).
- **Differentiate core taxa (circled) and discriminant taxa (solid points) in the legend**.

**Steps for Modifications:**

1. Extract taxonomy strings for core taxa (from `d70`).
2. Merge taxonomy info into the `data` dataframe.
3. Use `geom_text_repel()` from `ggrepel` to label core taxa without overlapping.

```r
# Load libraries
library(readxl)
library(dplyr)
library(ggplot2)
library(ggrepel)  # For non-overlapping text labels

# Read in main dataset
data <- read_excel("Maaslin_allresults_COMBINED.xlsx", sheet = "SF v TP") %>% as.data.frame()

# Read in discriminant taxa
discriminant_otus <- read_excel("Maaslin_discrim_results.xlsx", sheet = "SF5_discrim") %>%
  select(OTUID) %>%
  pull()

# Read in core taxa from d70 (core across all sites and depths)
core_taxa_data <- read_excel("CoreTaxa_70_16S_ASCC.xlsx", sheet = "d70") %>%
  select(OTUID, asv_string_silva)  # Keep taxonomy info

# Extract OTUIDs of core taxa
core_taxa <- core_taxa_data$OTUID

# Ensure matching format
data$OTUID <- trimws(as.character(data$OTUID))
discriminant_otus <- trimws(as.character(discriminant_otus))
core_taxa <- trimws(as.character(core_taxa))

# Merge taxonomy info into data
data <- data %>%
  left_join(core_taxa_data, by = "OTUID") %>%  # Add taxonomy string
  mutate(log_q = -log10(qval),
         Enrichment = case_when(
           qval >= 0.05 ~ "Non-Significant",
           coef > 0 ~ "SF Enriched",
           coef < 0 ~ "TP Enriched"
         ),
         Highlight = ifelse(qval < 0.05 & OTUID %in% discriminant_otus, "Highlighted", "Regular"),
         Core_Taxa = ifelse(OTUID %in% core_taxa, "Core", "Non-Core"))  # Mark core taxa

# Base volcano plot
plot <- ggplot(data, aes(x = coef, y = log_q, color = Enrichment)) +
  geom_point(size = 1.5, alpha = 0.9) +

  # Custom color scheme
  scale_color_manual(values = c(
    "SF Enriched" = "#CC79A7",
    "TP Enriched" = "#56B4E9",
    "Non-Significant" = "gray"
  )) +

  # Overlay highlighted discriminant taxa (black solid points)
  geom_point(data = data %>% filter(Highlight == "Highlighted"),
             aes(x = coef, y = log_q), color = "black", size = 1.5, alpha = .2, shape = 16) +

  # Overlay core taxa (circled points)
  geom_point(data = data %>% filter(Core_Taxa == "Core"),
             aes(x = coef, y = log_q), shape = 21, fill = NA, size = 2, stroke = 1, color = "black") +

  # Label core taxa with taxonomy strings
  geom_text_repel(data = data %>% filter(Core_Taxa == "Core"),
                  aes(x = coef, y = log_q, label = asv_string_silva),
                  size = 3, max.overlaps = 15, box.padding = 0.5, segment.color = "gray") +

  # Threshold lines
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "gray") +

  # Labels & theme
  labs(x = "Coefficient", y = "-log10(q-value)", color = "Enrichment") +
  theme_minimal() +
  theme(legend.position = "top")

# Show plot
print(plot)
```

</aside>

<aside>

## **Modifying the Volcano Plot to Label Core Taxa with Family-Level Taxonomy (no legend, using specific site core taxa)**

To label **core taxa** using the `"Family"` column instead of `"asv_string_silva"`, follow these steps.

---

### **1. Modify the Dataset to Add Family-Level Labels for Core Taxa**

***this is using the specific sheets for all taxa core for SF, and all taxa core to TP for one plots***

***We should also do this with the “d70” sheet***

```r
# Read core taxa and select the OTUID and Family column
core_taxa_sf <- read_excel("CoreTaxa_70_16S_ASCC.xlsx", sheet = "d70_SF") %>%
  select(OTUID, Family)

core_taxa_tp <- read_excel("CoreTaxa_70_16S_ASCC.xlsx", sheet = "d70_TP") %>%
  select(OTUID, Family)

# Merge taxonomy information into the main dataset
data <- data %>%
  left_join(bind_rows(core_taxa_sf, core_taxa_tp), by = "OTUID") %>%  # Add Family column
  mutate(Family_Label = case_when(
    Core_Taxa != "Non-Core" ~ Family,  # Assign Family names only to core taxa
    TRUE ~ NA_character_  # Keep all other points unlabeled
  ))

```

---

### **2. Modify the Volcano Plot to Label Core Taxa with Family Names**

```r
# Base volcano plot
plot <- ggplot(data, aes(x = coef, y = log_q, color = Enrichment)) +
  geom_point(size = 1.5, alpha = 0.9) +

  # Custom color scheme for enrichment
  scale_color_manual(values = c(
    "SF Enriched" = "#CC79A7",
    "TP Enriched" = "#56B4E9",
    "Non-Significant" = "gray"
  )) +

  # Overlay core taxa (circled and labeled)
  geom_point(data = data %>% filter(Core_Taxa != "Non-Core"),
             aes(x = coef, y = log_q), shape = 21, fill = NA, size = 2, stroke = 1, color = "black") +  # Black-outlined circles for core taxa

  # Add Family labels for core taxa
  geom_text_repel(data = data %>% filter(!is.na(Family_Label)),
                  aes(x = coef, y = log_q, label = Family_Label),
                  size = 3, color = "black", fontface = "italic", box.padding = 0.5, segment.color = "gray") +

  # Threshold lines
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "gray") +

  # Labels & theme
  labs(x = "Coefficient", y = "-log10(q-value)", color = "Enrichment") +
  theme_minimal() +
  theme(legend.position = "top")

# Show plot
print(plot)

```

---

## **Explanation of Changes**

- **Now labels core taxa with `"Family"` taxonomy instead of full taxonomy strings.**
- **Uses `geom_text_repel()`** to prevent label overlap and ensure clarity.
- **Core taxa remain circled (black-outlined, unfilled points).**
- **Only core taxa receive Family labels; all other taxa remain unlabeled.**

---

## **Final Thoughts**

This update:

- **Ensures core taxa are labeled with their Family names.**
- **Keeps the volcano plot clean while providing key taxonomic insights.**
- **Prevents text overlap for better readability.**

</aside>

<aside>

## **Adding a Legend for Core and Discriminant Taxa in the Volcano Plot**

To clearly **differentiate core taxa (circled) and discriminant taxa (solid points)** in the legend, update your plot with **custom shapes in the legend** while maintaining the existing plot structure.

---

### **1. Modify the Plot to Add a Legend Entry for Core and Discriminant Taxa**

```r
# Base volcano plot
plot <- ggplot(data, aes(x = coef, y = log_q, color = Enrichment)) +
  geom_point(size = 1.5, alpha = 0.9) +

  # Custom color scheme for enrichment
  scale_color_manual(values = c(
    "SF Enriched" = "#CC79A7",
    "TP Enriched" = "#56B4E9",
    "Non-Significant" = "gray"
  )) +

  # Overlay discriminant taxa (solid red points)
  geom_point(data = data %>% filter(Highlight == "Highlighted"),
             aes(x = coef, y = log_q, shape = "Discriminant Taxa"), size = 3, color = "#AB3329", alpha = 0.8) +

  # Overlay core taxa (circled points)
  geom_point(data = data %>% filter(Core_Taxa != "Non-Core"),
             aes(x = coef, y = log_q, shape = "Core Taxa"), fill = NA, size = 2, stroke = 1, color = "black") +

  # Define shape legend
  scale_shape_manual(values = c("Discriminant Taxa" = 16, "Core Taxa" = 21),
                     name = "Taxa Type") +

  # Threshold lines
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "gray") +

  # Labels & theme
  labs(x = "Coefficient", y = "-log10(q-value)", color = "Enrichment") +
  theme_minimal() +
  theme(legend.position = "top")

# Show plot
print(plot)

```

---

### **Explanation:**

- **Adds a custom legend** that labels **Discriminant Taxa (solid points)** and **Core Taxa (circled points)**.
- **Uses `scale_shape_manual()`** to define:
  - **Solid points (`shape = 16`) for discriminant taxa.**
  - **Circled points (`shape = 21`) for core taxa.**
- **Ensures clear legend labeling under "Taxa Type".**

---

## **Final Thoughts**

This update:

- **Adds a legend to distinguish between discriminant (solid) and core (circled) taxa.**
- **Ensures clarity in interpretation without adding extra labels to the plot.**
- **Keeps the visualization clean while maintaining ecological significance.**

</aside>

<aside>

## **Adding Labels for Highlighted Discriminant Points in the Volcano Plot (wouldn’t always want to do this)**

To label **highlighted discriminant points** (e.g., taxa that are significant in SF vs. TP and enriched at 0-5 cm) while keeping **circled points for core taxa**, follow these steps:

### **1. Modify the Dataset to Tag the Highlighted Discriminant Points**

```r
# Define a new category for special highlighted discriminant points
data <- data %>%
  mutate(Highlight_Label = case_when(
    Highlight == "Highlighted" ~ "Enriched in all sites at 0-5 cm",
    TRUE ~ NA_character_  # Keep all other points unlabeled
  ))

```

---

**Explanation:**

- Adds a **new column** (`Highlight_Label`) that assigns a label to **discriminant taxa**.
- All other points remain **unlabeled**.

---

### **2. Modify the Volcano Plot to Include Labels for Highlighted Discriminant Points**

```r
# Base volcano plot
plot <- ggplot(data, aes(x = coef, y = log_q, color = Enrichment)) +
  geom_point(size = 1.5, alpha = 0.9) +

  # Custom color scheme for enrichment
  scale_color_manual(values = c(
    "SF Enriched" = "#CC79A7",
    "TP Enriched" = "#56B4E9",
    "Non-Significant" = "gray"
  )) +

  # Overlay highlighted discriminant taxa (solid red points)
  geom_point(data = data %>% filter(Highlight == "Highlighted"),
             aes(x = coef, y = log_q), shape = 16, size = 3, color = "#AB3329", alpha = 0.8) +  # Solid red points for discriminant taxa

  # Overlay core taxa (circled points)
  geom_point(data = data %>% filter(Core_Taxa != "Non-Core"),
             aes(x = coef, y = log_q), shape = 21, fill = NA, size = 2, stroke = 1, color = "black") +  # Black-outlined circles for core taxa

  # Add labels for highlighted discriminant points
  geom_text_repel(data = data %>% filter(!is.na(Highlight_Label)),
                  aes(x = coef, y = log_q, label = Highlight_Label),
                  size = 4, color = "#AB3329", fontface = "bold", box.padding = 0.5, segment.color = "#AB3329") +

  # Threshold lines
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "gray") +

  # Labels & theme
  labs(x = "Coefficient", y = "-log10(q-value)", color = "Enrichment") +
  theme_minimal() +
  theme(legend.position = "top")

# Show plot
print(plot)

```

**Explanation:**

- **Discriminant taxa (highlighted points) are now solid red (`#AB3329`).**
- **Core taxa are circled (black-outlined, unfilled points).**
- **Labels the discriminant taxa** using `geom_text_repel()`, keeping labels non-overlapping.

---

## **Final Thoughts**

This update:

- **Highlights discriminant taxa in solid red (`#AB3329`).**
- **Keeps core taxa circled but unfilled.**
- **Adds dynamic labels for highlighted discriminant points** while preventing overlap.

This ensures clear distinction between:

- **Discriminant taxa driving differences (solid red).**
- **Core taxa found in multiple sites (circled black).**
- **Labeled significant points for better ecological interpretation.**

</aside>

<aside>

## **How to Incorporate This into Your Main Code**

To integrate **site-specific core taxa (d70_SF & d70_TP) and depth-specific core taxa (d70_shallow, d70_deep, d70_OM)** into your main R code, follow these steps:

---

### **1. Read Site-Specific Core Taxa (d70_SF & d70_TP)**

```r
# Read SF and TP core taxa separately
core_taxa_sf <- read_excel("CoreTaxa_70_16S_ASCC.xlsx", sheet = "d70_SF") %>%
  select(OTUID, asv_string_silva)

core_taxa_tp <- read_excel("CoreTaxa_70_16S_ASCC.xlsx", sheet = "d70_TP") %>%
  select(OTUID, asv_string_silva)

```

**Explanation:**

- Reads **core taxa for SF and TP separately**.
- Can be used to highlight **site-specific** core taxa in your volcano plot.

---

### **2. Combine Site-Specific Core Taxa (SF & TP)**

```r
# Merge SF and TP core taxa into one dataset
core_taxa_combined <- bind_rows(core_taxa_sf, core_taxa_tp) %>%
  distinct(OTUID, .keep_all = TRUE)

```

**Explanation:**

- Combines **both SF and TP core taxa** into one dataset.
- Ensures **each OTUID appears only once** to avoid duplication.

---

### **3. Read Depth-Specific Core Taxa (d70_shallow, d70_deep, d70_OM)**

```r
core_taxa_shallow <- read_excel("CoreTaxa_70_16S_ASCC.xlsx", sheet = "d70_shallow") %>%
  select(OTUID)

core_taxa_deep <- read_excel("CoreTaxa_70_16S_ASCC.xlsx", sheet = "d70_deep") %>%
  select(OTUID)

core_taxa_om <- read_excel("CoreTaxa_70_16S_ASCC.xlsx", sheet = "d70_OM") %>%
  select(OTUID)

```

**Explanation:**

- Reads **core taxa specific to different soil depths** (shallow, deep, and organic matter layers).
- This allows for **depth-specific visualization** in the volcano plot.

---

### **4. Incorporate Core Taxa into Your Main Dataset**

Modify your **main dataset** to flag site-specific and depth-specific core taxa:

```r
# Convert OTUIDs to character for consistency
data$OTUID <- trimws(as.character(data$OTUID))
core_taxa_sf$OTUID <- trimws(as.character(core_taxa_sf$OTUID))
core_taxa_tp$OTUID <- trimws(as.character(core_taxa_tp$OTUID))
core_taxa_combined$OTUID <- trimws(as.character(core_taxa_combined$OTUID))
core_taxa_shallow$OTUID <- trimws(as.character(core_taxa_shallow$OTUID))
core_taxa_deep$OTUID <- trimws(as.character(core_taxa_deep$OTUID))
core_taxa_om$OTUID <- trimws(as.character(core_taxa_om$OTUID))

# Add core taxa classifications to the dataset
data <- data %>%
  mutate(Core_Taxa = case_when(
    OTUID %in% core_taxa_sf$OTUID ~ "SF Core",
    OTUID %in% core_taxa_tp$OTUID ~ "TP Core",
    OTUID %in% core_taxa_shallow$OTUID ~ "Shallow Core",
    OTUID %in% core_taxa_deep$OTUID ~ "Deep Core",
    OTUID %in% core_taxa_om$OTUID ~ "OM Core",
    TRUE ~ "Non-Core"
  ))

```

**Explanation:**

- Ensures that **OTUIDs are formatted consistently** across datasets.
- Uses `mutate()` to **assign core taxonomy classifications**:
  - **"SF Core"** for taxa present in d70_SF.
  - **"TP Core"** for taxa present in d70_TP.
  - **"Shallow Core"**, **"Deep Core"**, and **"OM Core"** based on depth-specific sheets.
  - **"Non-Core"** for taxa that are not in any core dataset.

---

### **5. Modify Volcano Plot to Differentiate Core Taxa**

```r
# Base volcano plot
plot <- ggplot(data, aes(x = coef, y = log_q, color = Enrichment)) +
  geom_point(size = 1.5, alpha = 0.9) +

  # Custom color scheme for enrichment
  scale_color_manual(values = c(
    "SF Enriched" = "#CC79A7",
    "TP Enriched" = "#56B4E9",
    "Non-Significant" = "gray"
  )) +

  # Overlay site-specific and depth-specific core taxa with different markers
  geom_point(data = data %>% filter(Core_Taxa == "SF Core"),
             aes(x = coef, y = log_q), shape = 21, fill = NA, size = 2, stroke = 1, color = "#CC79A7") +  # Pink circle for SF core

  geom_point(data = data %>% filter(Core_Taxa == "TP Core"),
             aes(x = coef, y = log_q), shape = 22, fill = NA, size = 2, stroke = 1, color = "#56B4E9") +  # Blue square for TP core

  geom_point(data = data %>% filter(Core_Taxa == "Shallow Core"),
             aes(x = coef, y = log_q), shape = 23, fill = NA, size = 2, stroke = 1, color = "#E69F00") +  # Orange diamond for Shallow core

  geom_point(data = data %>% filter(Core_Taxa == "Deep Core"),
             aes(x = coef, y = log_q), shape = 24, fill = NA, size = 2, stroke = 1, color = "#009E73") +  # Green triangle for Deep core

  geom_point(data = data %>% filter(Core_Taxa == "OM Core"),
             aes(x = coef, y = log_q), shape = 25, fill = NA, size = 2, stroke = 1, color = "#F0E442") +  # Yellow inverted triangle for OM core

  # Label core taxa with taxonomy strings
  geom_text_repel(data = data %>% filter(Core_Taxa != "Non-Core"),
                  aes(x = coef, y = log_q, label = asv_string_silva),
                  size = 3, max.overlaps = 15, box.padding = 0.5, segment.color = "gray") +

  # Threshold lines
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "gray") +

  # Labels & theme
  labs(x = "Coefficient", y = "-log10(q-value)", color = "Enrichment") +
  theme_minimal() +
  theme(legend.position = "top")

# Show plot
print(plot)

```

**Explanation:**

- **Assigns unique shapes and colors** to different core taxa categories:
  - **Pink circles** for **SF Core**.
  - **Blue squares** for **TP Core**.
  - **Orange diamonds** for **Shallow Core**.
  - **Green triangles** for **Deep Core**.
  - **Yellow inverted triangles** for **OM Core**.
- **Labels core taxa** using `geom_text_repel()` to prevent overlap.

---

## **Final Thoughts**

By incorporating **site-specific and depth-specific core taxa**, this version of the volcano plot:

- **Differentiates core taxa by site and depth**.
- **Highlights key taxa driving differences between SF and TP**.
- **Allows better ecological interpretation of microbial distribution**.

This ensures your analysis **captures both site and depth-specific microbial patterns** while maintaining **statistical rigor** in identifying enriched taxa.

</aside>
