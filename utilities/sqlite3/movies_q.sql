-- SQL for the Actors table

BEGIN TRANSACTION;
CREATE TABLE Actors(AId integer primary key autoincrement, Name text);
INSERT INTO Actors VALUES(1,'Philip Seymour Hofman');
INSERT INTO Actors VALUES(2,'Kate Shindle');
INSERT INTO Actors VALUES(3,'Kelci Stephenson');
INSERT INTO Actors VALUES(4,'Al Pacino');
INSERT INTO Actors VALUES(5,'Gabrielle Anwar');
INSERT INTO Actors VALUES(6,'Patricia Arquette');
INSERT INTO Actors VALUES(7,'Gabriel Byrne');
INSERT INTO Actors VALUES(8,'Max von Sydow');
INSERT INTO Actors VALUES(9,'Ellen Burstyn');
INSERT INTO Actors VALUES(10,'Jason Miller');
COMMIT;
-- SQL for the Movies table

BEGIN TRANSACTION;
CREATE TABLE Movies(MId integer primary key autoincrement, Title text);
INSERT INTO Movies VALUES(1,'Capote');
INSERT INTO Movies VALUES(2,'Scent of a woman');
INSERT INTO Movies VALUES(3,'Stigmata');
INSERT INTO Movies VALUES(4,'Exorcist');
INSERT INTO Movies VALUES(5,'Hamsun');
COMMIT;
-- SQL for the ActorsMovies table

BEGIN TRANSACTION;
CREATE TABLE ActorsMovies(Id integer primary key autoincrement, 
AId integer, MId integer);
INSERT INTO ActorsMovies VALUES(1,1,1);
INSERT INTO ActorsMovies VALUES(2,2,1);
INSERT INTO ActorsMovies VALUES(3,3,1);
INSERT INTO ActorsMovies VALUES(4,4,2);
INSERT INTO ActorsMovies VALUES(5,5,2);
INSERT INTO ActorsMovies VALUES(6,6,3);
INSERT INTO ActorsMovies VALUES(7,7,3);
INSERT INTO ActorsMovies VALUES(8,8,4);
INSERT INTO ActorsMovies VALUES(9,9,4);
INSERT INTO ActorsMovies VALUES(10,10,4);
INSERT INTO ActorsMovies VALUES(11,8,5);
COMMIT;
