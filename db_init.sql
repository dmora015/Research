DROP TABLE Contains;
DROP TABLE Loci;
DROp TABLE Isolates;

CREATE TABLE Isolates	(id INTEGER UNIQUE NOT NULL,
			year INTEGER NOT NULL,
			name CHAR(5) NOT NULL,
			plant INTEGER NOT NULL,
			nodule INTEGER NOT NULL,
			site CHAR(30) NOT NULL,
			latitude REAL,
			longitude REAL,
			PRIMARY KEY(id),
			UNIQUE(year, name, plant, nodule));

CREATE TABLE Loci	(type ENUM('K','G','I','R','D','N','Z','L') NOT NULL,
			num INTEGER NOT NULL,
			bases TEXT NOT NULL,
			PRIMARY KEY(type, num));

CREATE TABLE Contains	(id INTEGER NOT NULL,
			type ENUM('K','G','I','R','D','N','Z','L') NOT NULL,
			num INTEGER NOT NULL,
			PRIMARY KEY(id, type),
			FOREIGN KEY(id) REFERENCES Isolates(id),
			FOREIGN KEY(type, num) REFERENCES Loci(type, num));