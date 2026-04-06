# Neo4j & UniProt: Exploring Human Protein Interactions

A hands-on tutorial for building and querying a human protein-protein interaction (PPI) graph database using Neo4j and UniProt data.

Graph databases excel at revealing hidden relationships in complex networks. Neo4j has been used across many fields — most famously in investigative journalism, where it played a central role in analyzing the **Panama Papers**: 11.5 million leaked financial documents exposing offshore tax havens. Journalists used Neo4j to map the connections between shell companies, intermediaries, and individuals that would have been impossible to detect in a traditional database. Read more: [Analyzing the Panama Papers with Neo4j](https://neo4j.com/blog/cypher-and-gql/analyzing-panama-papers-neo4j/)

The same graph thinking applies to biology — proteins interact in networks, diseases propagate through pathways, and hidden connections emerge when you query the right relationships.

---

## Table of Contents

1. [Why a Graph Database?](#1-why-a-graph-database)
2. [Installing Neo4j](#2-installing-neo4j)
3. [Downloading UniProt Data](#3-downloading-uniprot-data)
4. [Data Ingestion](#4-data-ingestion)
5. [Graph Model](#5-graph-model)
6. [Basic Queries](#6-basic-queries)
7. [Insulin Signaling Pathway](#7-insulin-signaling-pathway)
8. [Disease Queries](#8-disease-queries)

---

## 1. Why a Graph Database?

Protein-protein interaction data is inherently a network — proteins are nodes, interactions are edges. While you _could_ store this in a relational database, it quickly becomes painful:

| Question                           | Relational DB               | Neo4j                             |
| ---------------------------------- | --------------------------- | --------------------------------- |
| Who does protein A interact with?  | `JOIN` on interaction table | `MATCH (a)-[:INTERACTS_WITH]-(b)` |
| Find all proteins within 3 hops    | Recursive CTEs, slow        | `[:INTERACTS_WITH*..3]`, native   |
| Shortest path between two proteins | Complex, expensive          | `shortestPath()`, built-in        |
| Add a new relationship type        | Schema migration            | Just add a new relationship       |

The core issue with relational databases is that **path queries require self-joins** — each additional hop multiplies the query complexity. A 5-hop path query on a 20,000-protein dataset becomes impractical in SQL but is trivial in Cypher.

Graph databases also match how biologists think: _"What is between TLR4 and NF-κB?"_ is naturally a path query, not a table scan.

---

## 2. Installing Neo4j

See [INSTALL.md](INSTALL.md) for full installation instructions (Desktop, Server, Docker) and index creation.

---

## 3. Downloading UniProt Data

We use the **UniProt reviewed (Swiss-Prot)** database filtered for _Homo sapiens_ (taxon 9606).

### Steps

1. Go to [https://www.uniprot.org](https://www.uniprot.org)
2. Search for: `reviewed:true AND organism_id:9606`
3. Click **Download**
4. Choose format: **TSV**
5. Select the following columns (`Entry` is selected by default):

| Column                               | Description                                          |
| ------------------------------------ | ---------------------------------------------------- |
| `Entry`                              | UniProt accession (unique ID, e.g. `P04637`)         |
| `Entry Name`                         | Human-readable ID (e.g. `P53_HUMAN`)                 |
| `Protein names`                      | Full protein name                                    |
| `Gene Names`                         | Associated gene names                                |
| `Interacts with`                     | Interacting protein accessions (semicolon-separated) |
| `Involvement in disease`             | Associated diseases                                  |
| `Gene Ontology (biological process)` | GO biological process terms                          |
| `Gene Ontology (molecular function)` | GO molecular function terms                          |
| `Gene Ontology (cellular component)` | GO cellular component terms                          |

6. Download and save as `uniprotkb_human.tsv.gz`
7. Uncompress the file `gunzip uniprotkb_human.tsv.gz`
8. Check number of lines `wc -l uniprotkb_human.tsv` => 20432 uniprotkb_human.tsv

---

## 4. Data Ingestion

neo4j start

```cypher
// [PLACEHOLDER: run ingestion script here]
// The script loads uniprotkb_human.tsv, creates (:Protein) nodes
// and (:INTERACTS_WITH) relationships between them.
```

```bash
cat script.cql | cypher-shell -a bolt://localhost:7687 -u neo4j
```

The ingestion script handles:

- Creating one `(:Protein)` node per row using `Entry` as the unique ID
- Splitting the `Interacts with` column on `"; "` and creating `[:INTERACTS_WITH]` relationships
- Skipping null or empty interaction fields

---

## 5. Graph Model

```
(:Protein {
    entryID:      "P04637",        // UniProt accession
    entryNameID:  "P53_HUMAN",     // Entry name
    proteinName:  "Cellular tumor antigen p53",
    geneName:     "TP53",
    disease:      "Li-Fraumeni syndrome; ..."
})-[:INTERACTS_WITH]->(:Protein)
```

Relationships are directional as loaded from UniProt, but are queried without direction (undirected) unless otherwise specified.

---

## 6. Basic Queries

### 5.1 Direct interaction partners

Who does a given protein directly interact with?

```cypher
MATCH (p:Protein {entryNameID: 'INSR_HUMAN'})-[:INTERACTS_WITH]-(partner:Protein)
RETURN partner.entryNameID, partner.proteinName
```

### 5.2 Most connected proteins (hubs)

Which proteins have the most interactions in the dataset?

```cypher
MATCH (p:Protein)-[:INTERACTS_WITH]-(partner:Protein)
RETURN p.entryNameID, count(partner) AS degree
ORDER BY degree DESC
LIMIT 10
```

> Hub proteins are often chaperones, ubiquitin ligases, or scaffolding proteins — biologically important but can create noise in pathway queries.

### 5.3 Shared interaction partners

Do two proteins share common interactors?

```cypher
MATCH (p1:Protein {entryNameID: 'INSR_HUMAN'})-[:INTERACTS_WITH]-(common:Protein)-[:INTERACTS_WITH]-(p2:Protein {entryNameID: 'IRS1_HUMAN'})
RETURN common.entryNameID AS sharedPartner
```

### 5.4 Shortest path between two proteins

```cypher
MATCH path = shortestPath(
  (p1:Protein {entryNameID: 'INSR_HUMAN'})-[:INTERACTS_WITH*..5]-(p2:Protein {entryNameID: 'AKT1_HUMAN'})
)
RETURN [n IN nodes(path) | n.entryNameID] AS pathway, length(path) AS hops
```

### 5.5 All shortest paths between two proteins

```cypher
MATCH path = allShortestPaths(
  (p1:Protein {entryNameID: 'INSR_HUMAN'})-[:INTERACTS_WITH*..5]-(p2:Protein {entryNameID: 'AKT1_HUMAN'})
)
RETURN [n IN nodes(path) | n.entryNameID] AS pathway, length(path) AS hops
ORDER BY hops
```

### 5.6 Explore a local subgraph

Visualize the interaction neighborhood up to 3 hops from a protein — best run in Neo4j Browser for graph visualization.

```cypher
MATCH path = (p:Protein {entryNameID: 'INSR_HUMAN'})-[:INTERACTS_WITH*..3]-(partner:Protein)
RETURN path
```

---

## 7. Insulin Signaling Pathway

The insulin signaling pathway is a textbook example well suited for graph exploration.

### Key proteins

| Protein          | entryNameID    | Role                           |
| ---------------- | -------------- | ------------------------------ |
| Insulin receptor | `INSR_HUMAN`   | Receptor tyrosine kinase       |
| IRS1             | `IRS1_HUMAN`   | Docking/adaptor protein        |
| PI3K regulatory  | `PIK3R1_HUMAN` | Lipid kinase complex           |
| PI3K catalytic   | `PIK3CA_HUMAN` | Lipid kinase complex           |
| AKT1             | `AKT1_HUMAN`   | Serine/threonine kinase        |
| mTOR             | `MTOR_HUMAN`   | Growth & metabolism regulator  |
| GLUT4            | `SLC2A4_HUMAN` | Glucose transporter            |
| GRB2             | `GRB2_HUMAN`   | Adaptor, links to MAPK cascade |

### 6.1 What does the insulin receptor interact with?

```cypher
MATCH (p:Protein {entryNameID: 'INSR_HUMAN'})-[:INTERACTS_WITH]-(partner:Protein)
RETURN partner.entryNameID, partner.proteinName
```

### 6.2 What does IRS1 interact with?

```cypher
MATCH (p:Protein {entryNameID: 'IRS1_HUMAN'})-[:INTERACTS_WITH]-(partner:Protein)
RETURN partner.entryNameID, partner.proteinName
```

### 6.3 Shared partners between INSR and IRS1

```cypher
MATCH (p1:Protein {entryNameID: 'INSR_HUMAN'})-[:INTERACTS_WITH]-(common:Protein)-[:INTERACTS_WITH]-(p2:Protein {entryNameID: 'IRS1_HUMAN'})
RETURN common.entryNameID AS sharedPartner
```

### 6.4 Shortest path from INSR to AKT1

```cypher
MATCH path = shortestPath(
  (p1:Protein {entryNameID: 'INSR_HUMAN'})-[:INTERACTS_WITH*..5]-(p2:Protein {entryNameID: 'AKT1_HUMAN'})
)
RETURN [n IN nodes(path) | n.entryNameID] AS pathway, length(path) AS hops
```

### 6.5 All shortest paths from INSR to AKT1

```cypher
MATCH path = allShortestPaths(
  (p1:Protein {entryNameID: 'INSR_HUMAN'})-[:INTERACTS_WITH*..5]-(p2:Protein {entryNameID: 'AKT1_HUMAN'})
)
RETURN [n IN nodes(path) | n.entryNameID] AS pathway, length(path) AS hops
ORDER BY hops
```

### 6.6 Most connected proteins in the insulin subgraph

Which proteins are the most central within 3 hops of INSR?

```cypher
MATCH (p:Protein {entryNameID: 'INSR_HUMAN'})-[:INTERACTS_WITH*..3]-(partner:Protein)
WITH partner, count(*) AS connections
RETURN partner.entryNameID, connections
ORDER BY connections DESC
LIMIT 10
```

> Hub proteins like PIK3R1, AKT1, and GRB2 naturally emerge here — they are genuine pathway members but also interact with many other proteins.

---

## 8. Disease Queries

UniProt includes disease association annotations for many proteins. These enable clinically relevant graph queries.

### 7.1 Proteins associated with a disease

Find all proteins linked to Type 2 diabetes:

```cypher
MATCH (p:Protein)
WHERE p.disease CONTAINS 'diabetes'
RETURN p.entryNameID, p.proteinName, p.disease
```

### 7.2 Cancer-associated proteins that interact with p53

```cypher
MATCH (p53:Protein {entryNameID: 'P53_HUMAN'})-[:INTERACTS_WITH]-(partner:Protein)
WHERE partner.disease CONTAINS 'cancer'
   OR partner.disease CONTAINS 'carcinoma'
   OR partner.disease CONTAINS 'tumor'
RETURN partner.entryNameID, partner.proteinName, partner.disease
```

### 7.3 Shortest path between two disease-associated proteins

Path between the insulin receptor (Type 2 diabetes) and BRCA1 (breast cancer):

```cypher
MATCH path = shortestPath(
  (p1:Protein {entryNameID: 'INSR_HUMAN'})-[:INTERACTS_WITH*..6]-(p2:Protein {entryNameID: 'BRCA1_HUMAN'})
)
RETURN [n IN nodes(path) | n.entryNameID] AS pathway, length(path) AS hops
```

### 7.4 Proteins involved in multiple diseases

```cypher
MATCH (p:Protein)
WHERE p.disease IS NOT NULL AND p.disease <> ''
RETURN p.entryNameID, p.proteinName, size(split(p.disease, ';')) AS diseaseCount, p.disease
ORDER BY diseaseCount DESC
LIMIT 10
```

### 7.5 Interaction network around a disease gene

Explore the neighborhood of BRCA1 — useful for identifying candidate drug targets:

```cypher
MATCH path = (p:Protein {entryNameID: 'BRCA1_HUMAN'})-[:INTERACTS_WITH*..2]-(partner:Protein)
RETURN path
```

### 7.6 Disease proteins that are hubs

Which disease-associated proteins have the most interactions? Potential drug targets.

```cypher
MATCH (p:Protein)-[:INTERACTS_WITH]-(partner:Protein)
WHERE p.disease IS NOT NULL AND p.disease <> ''
RETURN p.entryNameID, p.proteinName, count(partner) AS degree, p.disease
ORDER BY degree DESC
LIMIT 10
```

### 7.7 Do two disease proteins share interaction partners?

Common partners between TP53 (cancer) and INSR (diabetes) may represent crosstalk proteins:

```cypher
MATCH (p1:Protein {entryNameID: 'P53_HUMAN'})-[:INTERACTS_WITH]-(common:Protein)-[:INTERACTS_WITH]-(p2:Protein {entryNameID: 'INSR_HUMAN'})
RETURN common.entryNameID AS sharedPartner, common.proteinName
```

---

## Key Concepts Covered

| Concept                | Cypher Feature                                   |
| ---------------------- | ------------------------------------------------ |
| Single node lookup     | `MATCH (p:Protein {entryNameID: '...'})`         |
| Direct neighbors       | `-[:INTERACTS_WITH]-`                            |
| Aggregation & ranking  | `count()`, `ORDER BY`, `LIMIT`                   |
| 2-hop patterns         | `-[:INTERACTS_WITH]-(common)-[:INTERACTS_WITH]-` |
| Shortest path          | `shortestPath()` with `*..n` bound               |
| All shortest paths     | `allShortestPaths()` with `*..n` bound           |
| Subgraph visualization | `RETURN path` in Neo4j Browser                   |
| Text filtering         | `CONTAINS`, `IS NOT NULL`                        |
| List operations        | `split()`, `size()`                              |

---

## Further Reading

- UniProt documentation: [https://www.uniprot.org/help](https://www.uniprot.org/help)
- Cypher query language: [https://neo4j.com/docs/cypher-manual](https://neo4j.com/docs/cypher-manual)
- KEGG TLR signaling pathway: `hsa04620`
- Reactome insulin signaling: [https://reactome.org](https://reactome.org)
