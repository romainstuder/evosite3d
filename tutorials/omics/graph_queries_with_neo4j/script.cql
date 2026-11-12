// Delete all indices and all nodes / edges
DROP INDEX protein_index_name IF EXISTS;
DROP INDEX disease_index_name IF EXISTS;
DROP INDEX go_bp_index_name IF EXISTS;
DROP INDEX go_mf_index_name IF EXISTS;
DROP INDEX go_cc_index_name IF EXISTS;

MATCH (n)
OPTIONAL MATCH (n)-[r]-()
DELETE n, r;

// Create index for fast lookups
CREATE INDEX protein_index_name FOR (p:Protein) ON (p.entryID);
CREATE INDEX disease_index_name FOR (d:Disease) ON (d.id);
CREATE INDEX go_bp_index_name FOR (bp:BiologicalProcess) ON (bp.id);
CREATE INDEX go_mf_index_name FOR (mf:MolecularFunction) ON (mf.id);
CREATE INDEX go_cc_index_name FOR (cc:CellularComponent) ON (cc.id);

// Import data and create connections
WITH "file:///Users/romainstuder/Github/evosite3d/tutorials/omics/graph_queries_with_neo4j/uniprotkb_human.tsv" AS file_path
LOAD CSV WITH HEADERS FROM file_path AS row
FIELDTERMINATOR '\t'
MERGE (p1:Protein {entryID:row.Entry})
ON CREATE SET p1.entryNameID = row.`Entry Name`,
              p1.proteinName = row.`Protein names`,
              p1.geneName = row.`Gene Names`,
              p1.diseaseRaw = row.`Involvement in disease`,
              p1.goBiologicalProcessRaw = row.`Gene Ontology (biological process)`,
              p1.goMolecularFunctionRaw = row.`Gene Ontology (molecular function)`,
              p1.goCellularComponentRaw = row.`Gene Ontology (cellular component)`
WITH p1, split(row.`Interacts with`, ';') AS partners
UNWIND partners AS partnerId
MATCH (p2:Protein {entryID: trim(partnerId)})
MERGE (p1)-[:INTERACTS_WITH]->(p2);

// Create Disease nodes from protein disease annotations
// Format: "DISEASE: <name> [MIM:<id>]: <description>..."
WITH "file:///Users/romainstuder/Github/evosite3d/tutorials/omics/graph_queries_with_neo4j/uniprotkb_human.tsv" AS file_path
LOAD CSV WITH HEADERS FROM file_path AS row
FIELDTERMINATOR '\t'
WITH row
WHERE row.`Involvement in disease` IS NOT NULL AND row.`Involvement in disease` <> ''
MATCH (p:Protein {entryID: row.Entry})
WITH p, split(row.`Involvement in disease`, 'DISEASE: ') AS diseases
UNWIND diseases AS entry
WITH p, trim(entry) AS entry
WHERE entry <> '' AND entry CONTAINS '[MIM:'
WITH p,
     trim(split(entry, '[MIM:')[0]) AS diseaseName,
     split(split(entry, '[MIM:')[1], ']')[0] AS mimId
WHERE diseaseName <> '' AND mimId <> ''
MERGE (d:Disease {id: mimId})
ON CREATE SET d.name = diseaseName
MERGE (p)-[:ASSOCIATED_WITH]->(d);

// Create BiologicalProcess nodes from GO biological process annotations
WITH "file:///Users/romainstuder/Github/evosite3d/tutorials/omics/graph_queries_with_neo4j/uniprotkb_human.tsv" AS file_path
LOAD CSV WITH HEADERS FROM file_path AS row
FIELDTERMINATOR '\t'
WITH row
WHERE row.`Gene Ontology (biological process)` IS NOT NULL AND row.`Gene Ontology (biological process)` <> ''
MATCH (p:Protein {entryID: row.Entry})
WITH p, split(row.`Gene Ontology (biological process)`, ';') AS terms
UNWIND terms AS term
WITH p, trim(term) AS term
WHERE term <> ''
WITH p, trim(split(term, '[')[0]) AS goName, replace(split(term, '[')[1], ']', '') AS goId
MERGE (bp:BiologicalProcess {id: goId})
ON CREATE SET bp.name = goName
MERGE (p)-[:INVOLVED_IN]->(bp);

// Create MolecularFunction nodes from GO molecular function annotations
WITH "file:///Users/romainstuder/Github/evosite3d/tutorials/omics/graph_queries_with_neo4j/uniprotkb_human.tsv" AS file_path
LOAD CSV WITH HEADERS FROM file_path AS row
FIELDTERMINATOR '\t'
WITH row
WHERE row.`Gene Ontology (molecular function)` IS NOT NULL AND row.`Gene Ontology (molecular function)` <> ''
MATCH (p:Protein {entryID: row.Entry})
WITH p, split(row.`Gene Ontology (molecular function)`, ';') AS terms
UNWIND terms AS term
WITH p, trim(term) AS term
WHERE term <> ''
WITH p, trim(split(term, '[')[0]) AS goName, replace(split(term, '[')[1], ']', '') AS goId
MERGE (mf:MolecularFunction {id: goId})
ON CREATE SET mf.name = goName
MERGE (p)-[:HAS_FUNCTION]->(mf);

// Create CellularComponent nodes from GO cellular component annotations
WITH "file:///Users/romainstuder/Github/evosite3d/tutorials/omics/graph_queries_with_neo4j/uniprotkb_human.tsv" AS file_path
LOAD CSV WITH HEADERS FROM file_path AS row
FIELDTERMINATOR '\t'
WITH row
WHERE row.`Gene Ontology (cellular component)` IS NOT NULL AND row.`Gene Ontology (cellular component)` <> ''
MATCH (p:Protein {entryID: row.Entry})
WITH p, split(row.`Gene Ontology (cellular component)`, ';') AS terms
UNWIND terms AS term
WITH p, trim(term) AS term
WHERE term <> ''
WITH p, trim(split(term, '[')[0]) AS goName, replace(split(term, '[')[1], ']', '') AS goId
MERGE (cc:CellularComponent {id: goId})
ON CREATE SET cc.name = goName
MERGE (p)-[:IS_LOCATED_IN]->(cc);
