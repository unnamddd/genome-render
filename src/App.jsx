import React, { useEffect, useRef, useState } from "react";
import * as THREE from "three";
import { OrbitControls } from "three/addons/controls/OrbitControls.js";

const GenomeViewer = () => {
  const mountRef = useRef(null);
  const [loading, setLoading] = useState(true);
  const [genomeData, setGenomeData] = useState(null);
  const [selectedFeature, setSelectedFeature] = useState(null);
  const [metadata, setMetadata] = useState(null);
  const [viewMode, setViewMode] = useState('structure'); // 'genome', 'structure', or 'combined'
  const [structureData, setStructureData] = useState(null);

  useEffect(() => {
    let scene, camera, renderer, genomeObject, structureObject;
    let features = [];
    let annotationLabels = [];

    const setupControls = () => {
      console.log(THREE);
      const controls = new OrbitControls(camera, renderer.domElement);

      // Configure controls
      controls.enableDamping = true; // Smooth camera movements
      controls.dampingFactor = 0.05;

      // Set zoom limits
      controls.minDistance = 20; // Minimum zoom distance
      controls.maxDistance = 200; // Maximum zoom distance

      // Limit vertical rotation if desired
      controls.minPolar = Math.PI / 4; // Limit how high user can orbit
      controls.maxPolar = (Math.PI * 3) / 4; // Limit how low user can orbit

      // Optional: Disable vertical rotation entirely
      // controls.enableRotate = false;

      // Enable panning
      controls.enablePan = true;
      controls.panSpeed = 0.5; // Adjust pan sensitivity

      // Prevent zooming getting too close or far
      controls.minZoom = 0.5;
      controls.maxZoom = 2;

      return controls;
    };

    const parseGenBankData = (data) => {
      const features = [];
      const metadata = {
        accession: '',
        definition: '',
        organism: '',
        version: '',
        keywords: [],
        source: {},
        references: []
      };

      let currentFeature = null;
      let currentSection = '';

      const lines = data.split('\n');
      console.log(`Total lines to parse: ${lines.length}`);

      lines.forEach((line, index) => {
        // Feature parsing
        if (line.startsWith('     ')) {
          // New feature
          if (line.substring(5, 21).trim()) {
            if (currentFeature) {
              features.push(currentFeature);
            }

            const featureType = line.substring(5, 21).trim();
            const locationStr = line.substring(21).trim();
            console.log(`Parsing new feature: ${featureType} at location: ${locationStr}`);

            // Parse the location string
            let parsedLocation = parseLocation(locationStr);
            console.log('Parsed location result:', parsedLocation);

            currentFeature = {
              type: featureType,
              rawLocation: locationStr,  // Store raw location string for debugging
              location: parsedLocation,
              properties: {},
              coordinates: parsedLocation, // Directly use parsed location as coordinates
              strand: locationStr.includes('complement') ? '-' : '+'
            };
          }
          // Feature property
          else if (currentFeature && line.includes('/')) {
            const [key, ...valueParts] = line.trim().substring(1).split('=');
            const value = valueParts.join('=').replace(/"/g, '');
            currentFeature.properties[key] = value;
          }
        }
      });

      if (currentFeature) {
        features.push(currentFeature);
      }

      return { features, metadata };
    };

    const parseLocation = (location) => {
      console.log('Starting location parse:', location);

      // Function to extract numbers from a position string
      const extractPosition = (pos) => {
        const num = parseInt(pos.replace(/[^0-9]/g, ''));
        console.log(`Extracted position ${num} from ${pos}`);
        return num;
      };

      // Function to parse a single range or position
      const parseRange = (range) => {
        console.log('Parsing range:', range);
        if (range.includes('..')) {
          const [start, end] = range.split('..');
          return [extractPosition(start), extractPosition(end)];
        } else {
          const pos = extractPosition(range);
          return [pos, pos];
        }
      };

      try {
        // Remove any whitespace
        location = location.trim();

        // Store original location for debugging
        console.log('Processing location:', location);

        // Remove complement wrapper if present
        if (location.startsWith('complement(')) {
          location = location.match(/complement\((.*)\)/)[1];
          console.log('Removed complement wrapper:', location);
        }

        // Handle join operations
        if (location.startsWith('join(')) {
          location = location.match(/join\((.*)\)/)[1];
          console.log('Processing joined location:', location);

          // Split segments and parse each one
          const segments = location.split(',').map(s => s.trim());
          console.log('Split into segments:', segments);

          return segments.map(segment => parseRange(segment));
        }

        // Handle simple ranges or single positions
        return [parseRange(location)];

      } catch (error) {
        console.error('Error parsing location:', location, error);
        return [[0, 0]]; // Return a default value instead of empty array
      }
    };

    const createGenomeVisualization = (data) => {
      const genomeGroup = new THREE.Group();
      const features = [];

      // Configuration
      const genomeLength = data.metadata.length || 30000;
      const visualLength = 100; // Total visual length
      const scale = visualLength / genomeLength;

      // Create DNA double helix backbone
      const createDoubleHelix = () => {
        const helixGroup = new THREE.Group();
        const points1 = [];
        const points2 = [];

        // Parameters for double helix
        const radius = 1;  // Radius of helix
        const pitch = 5;   // Distance between turns
        const turns = 50;  // Number of turns
        const pointsPerTurn = 20;

        // Create points for both strands
        for (let i = 0; i <= turns * pointsPerTurn; i++) {
          const t = i / pointsPerTurn;
          const angle = t * Math.PI * 2;

          // First strand
          points1.push(new THREE.Vector3(
            t * pitch,
            Math.cos(angle) * radius,
            Math.sin(angle) * radius
          ));

          // Second strand (offset by 180 degrees)
          points2.push(new THREE.Vector3(
            t * pitch,
            Math.cos(angle + Math.PI) * radius,
            Math.sin(angle + Math.PI) * radius
          ));
        }

        // Create tubes for both strands
        const tubeGeometry1 = new THREE.TubeGeometry(
          new THREE.CatmullRomCurve3(points1),
          turns * pointsPerTurn,
          0.2,
          8,
          false
        );

        const tubeGeometry2 = new THREE.TubeGeometry(
          new THREE.CatmullRomCurve3(points2),
          turns * pointsPerTurn,
          0.2,
          8,
          false
        );

        const strandMaterial = new THREE.MeshPhongMaterial({
          color: 0xcccccc,
          transparent: true,
          opacity: 0.8,
          shininess: 30
        });

        const strand1 = new THREE.Mesh(tubeGeometry1, strandMaterial);
        const strand2 = new THREE.Mesh(tubeGeometry2, strandMaterial);

        helixGroup.add(strand1);
        helixGroup.add(strand2);

        return helixGroup;
      };

      // Create base DNA structure
      const dnaStructure = createDoubleHelix();
      genomeGroup.add(dnaStructure);

      // Function to map genome position to helix position
      const mapToHelix = (position) => {
        const normalizedPos = position / genomeLength;
        const turns = 50; // Should match the number in createDoubleHelix
        const pitch = 5;  // Should match the pitch in createDoubleHelix

        const t = normalizedPos * turns;
        const angle = t * Math.PI * 2;

        return {
          x: t * pitch,
          y: Math.cos(angle),
          z: Math.sin(angle)
        };
      };

      // Create features along the helix
      data.features.forEach((feature, index) => {
        if (!feature.coordinates || !feature.coordinates.length) {
          console.warn(`Feature ${index} has no coordinates:`, feature);
          return;
        }

        feature.coordinates.forEach(coordSet => {
          const [start, end] = coordSet;

          if (!start || !end || isNaN(start) || isNaN(end)) {
            console.warn(`Invalid coordinates for feature ${index}:`, coordSet);
            return;
          }

          try {
            // Get positions along helix
            const startPos = mapToHelix(start);
            const endPos = mapToHelix(end);

            // Create feature visualization
            let geometry, material;

            switch (feature.type.toLowerCase()) {
              case 'gene':
              case 'cds': {
                // Create curved feature following helix
                const points = [];
                const steps = 10;

                for (let i = 0; i <= steps; i++) {
                  const t = i / steps;
                  const pos = mapToHelix(start + (end - start) * t);
                  points.push(new THREE.Vector3(pos.x, pos.y * 1.5, pos.z * 1.5));
                }

                const curve = new THREE.CatmullRomCurve3(points);
                geometry = new THREE.TubeGeometry(curve, 8, 0.3, 8, false);
                material = new THREE.MeshPhongMaterial({
                  color: feature.type.toLowerCase() === 'gene' ? 0x00ff00 : 0xff0000,
                  transparent: true,
                  opacity: 0.8
                });
                break;
              }

              default: {
                // Create marker at feature position
                geometry = new THREE.SphereGeometry(0.4);
                material = new THREE.MeshPhongMaterial({
                  color: 0x0000ff,
                  transparent: true,
                  opacity: 0.8
                });
                const mesh = new THREE.Mesh(geometry, material);
                mesh.position.set(startPos.x, startPos.y * 1.5, startPos.z * 1.5);
                break;
              }
            }

            const featureMesh = new THREE.Mesh(geometry, material);
            featureMesh.userData = feature;
            features.push(featureMesh);
            genomeGroup.add(featureMesh);

          } catch (error) {
            console.error('Error creating feature visualization:', error);
          }
        });
      });

      genomeGroup.position.set(-visualLength / 2, 0, 0);
      scene.add(genomeGroup);
      return { genomeGroup };
    };


    const fetchGenomeData = async () => {
      try {
        const response = await fetch(
          "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=MW004169&rettype=gb&retmode=text"
        );
        const text = await response.text();

        console.log(
          "Raw GenBank data first 500 chars:",
          text.substring(0, 500)
        );

        // Parse GenBank format
        const parsedData = parseGenBankData(text);
        console.log("Parsed genome data:", {
          totalFeatures: parsedData.features.length,
          metadata: parsedData.metadata,
          sampleFeatures: parsedData.features.slice(0, 3), // First 3 features for inspection
        });

        setGenomeData(parsedData);
        return parsedData;
      } catch (error) {
        console.error("Error fetching genome data:", error);
        return null;
      }
    };

    const fetchPDBData = async (pdbId) => {
      try {
        // Fetch from PDB's REST API
        const response = await fetch(`https://files.rcsb.org/download/${pdbId}.pdb`);
        if (!response.ok) {
          throw new Error(`Failed to fetch PDB data: ${response.statusText}`);
        }
        const pdbData = await response.text();
        console.log(`Fetched PDB data for ${pdbId}, length: ${pdbData.length}`);
        return pdbData;
      } catch (error) {
        console.error('Error fetching PDB data:', error);
        return null;
      }
    };

    const parsePDBStructure = (pdbData) => {
      const atoms = [];
      const helices = [];
      const sheets = [];
      const models = [];

      const lines = pdbData.split('\n');
      let currentModel = null;

      lines.forEach(line => {
        const recordType = line.substring(0, 6).trim();

        switch (recordType) {
          case 'ATOM':
          case 'HETATM': {
            const atom = {
              serial: parseInt(line.substring(6, 11)),
              name: line.substring(12, 16).trim(),
              altLoc: line.substring(16, 17),
              resName: line.substring(17, 20).trim(),
              chainID: line.substring(21, 22),
              resSeq: parseInt(line.substring(22, 26)),
              iCode: line.substring(26, 27),
              x: parseFloat(line.substring(30, 38)),
              y: parseFloat(line.substring(38, 46)),
              z: parseFloat(line.substring(46, 54)),
              occupancy: parseFloat(line.substring(54, 60)),
              tempFactor: parseFloat(line.substring(60, 66)),
              element: line.substring(76, 78).trim(),
              charge: line.substring(78, 80).trim()
            };

            if (currentModel) {
              currentModel.atoms.push(atom);
            } else {
              atoms.push(atom);
            }
            break;
          }

          case 'HELIX': {
            helices.push({
              serNum: parseInt(line.substring(7, 10)),
              helixID: line.substring(11, 14).trim(),
              initResName: line.substring(15, 18).trim(),
              initChainID: line.substring(19, 20),
              initSeqNum: parseInt(line.substring(21, 25)),
              initICode: line.substring(25, 26),
              endResName: line.substring(27, 30).trim(),
              endChainID: line.substring(31, 32),
              endSeqNum: parseInt(line.substring(33, 37)),
              endICode: line.substring(37, 38),
              helixClass: parseInt(line.substring(38, 40)),
              length: parseInt(line.substring(71, 76))
            });
            break;
          }

          case 'SHEET': {
            sheets.push({
              strand: parseInt(line.substring(7, 10)),
              sheetID: line.substring(11, 14).trim(),
              numStrands: parseInt(line.substring(14, 16)),
              initResName: line.substring(17, 20).trim(),
              initChainID: line.substring(21, 22),
              initSeqNum: parseInt(line.substring(22, 26)),
              initICode: line.substring(26, 27),
              endResName: line.substring(28, 31).trim(),
              endChainID: line.substring(32, 33),
              endSeqNum: parseInt(line.substring(33, 37)),
              endICode: line.substring(37, 38)
            });
            break;
          }

          case 'MODEL': {
            currentModel = {
              serial: parseInt(line.substring(10, 14)),
              atoms: []
            };
            break;
          }

          case 'ENDMDL': {
            if (currentModel) {
              models.push(currentModel);
              currentModel = null;
            }
            break;
          }
        }
      });

      return {
        atoms,
        helices,
        sheets,
        models: models.length > 0 ? models : [{ serial: 1, atoms }]
      };
    };

    const createAtomMesh = (atom, offset) => {
      const radius = getAtomRadius(atom.element);
      const geometry = new THREE.SphereGeometry(radius, 16, 16);
      const material = new THREE.MeshPhongMaterial({
        color: getAtomColor(atom.element),
        shininess: 30,
        transparent: true,
        opacity: 0.8
      });

      const mesh = new THREE.Mesh(geometry, material);
      mesh.position.set(
        atom.x + offset.x,
        atom.y + offset.y,
        atom.z + offset.z
      );

      // Store atom data for interaction
      mesh.userData = {
        type: 'atom',
        element: atom.element,
        name: atom.name,
        resName: atom.resName,
        resSeq: atom.resSeq,
        chainID: atom.chainID
      };

      return mesh;
    };

    const createBond = (atom1, atom2, offset) => {
      const start = new THREE.Vector3(
        atom1.x + offset.x,
        atom1.y + offset.y,
        atom1.z + offset.z
      );
      const end = new THREE.Vector3(
        atom2.x + offset.x,
        atom2.y + offset.y,
        atom2.z + offset.z
      );

      const bondGeom = createCylinderBetweenPoints(start, end, 0.1);
      const bondMaterial = new THREE.MeshPhongMaterial({
        color: 0x808080,
        opacity: 0.5,
        transparent: true
      });

      return new THREE.Mesh(bondGeom, bondMaterial);
    };

    const createCylinderBetweenPoints = (start, end, radius) => {
      const direction = new THREE.Vector3().subVectors(end, start);
      const length = direction.length();

      const geometry = new THREE.CylinderGeometry(radius, radius, length, 8, 1);
      geometry.translate(0, length / 2, 0);
      geometry.rotateX(Math.PI / 2);

      const quaternion = new THREE.Quaternion();
      quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), direction.normalize());

      geometry.applyQuaternion(quaternion);
      geometry.translate(start.x, start.y, start.z);

      return geometry;
    };

    const getAtomRadius = (element) => {
      const radii = {
        'H': 0.1,
        'C': 0.2,
        'N': 0.2,
        'O': 0.2,
        'P': 0.25,
        'S': 0.25,
        'default': 0.2
      };
      return radii[element] || radii.default;
    };

    const getAtomColor = (element) => {
      const colors = {
        'H': 0xFFFFFF,
        'C': 0x808080,
        'N': 0x0000FF,
        'O': 0xFF0000,
        'P': 0xFFA500,
        'S': 0xFFFF00,
        'default': 0x008000
      };
      return colors[element] || colors.default;
    };

    const calculateDistance = (atom1, atom2) => {
      return Math.sqrt(
        Math.pow(atom2.x - atom1.x, 2) +
        Math.pow(atom2.y - atom1.y, 2) +
        Math.pow(atom2.z - atom1.z, 2)
      );
    };

    // Color generator for chains
    const ColorGenerator = function() {
      const colors = [
        0xFF0000, 0x00FF00, 0x0000FF,
        0xFFFF00, 0xFF00FF, 0x00FFFF,
        0xFF8000, 0x8000FF, 0x00FF80
      ];
      let index = 0;

      this.getNextColor = function() {
        const color = colors[index];
        index = (index + 1) % colors.length;
        return color;
      };
    };

    const createPDBVisualization = (structure, scene) => {
      console.log('Creating PDB visualization');
      const structureGroup = new THREE.Group();

      // Track min/max for centering
      let minX = Infinity, minY = Infinity, minZ = Infinity;
      let maxX = -Infinity, maxY = -Infinity, maxZ = -Infinity;

      // First pass: calculate exact bounds
      structure.models[0].atoms.forEach(atom => {
        minX = Math.min(minX, atom.x);
        minY = Math.min(minY, atom.y);
        minZ = Math.min(minZ, atom.z);
        maxX = Math.max(maxX, atom.x);
        maxY = Math.max(maxY, atom.y);
        maxZ = Math.max(maxZ, atom.z);
      });

      // Calculate exact center
      const centerOffset = {
        x: -(minX + maxX) / 2 - 100,
        y: -(minY + maxY) / 2,
        z: -(minZ + maxZ) / 2
      };

      console.log('Calculated center:', centerOffset);

      // Group atoms by chain and residue
      const chains = new Map();
      structure.models[0].atoms.forEach(atom => {
        if (!chains.has(atom.chainID)) {
          chains.set(atom.chainID, new Map());
        }
        const chain = chains.get(atom.chainID);

        if (!chain.has(atom.resSeq)) {
          chain.set(atom.resSeq, []);
        }
        chain.get(atom.resSeq).push(atom);
      });

      // Create chain groups with different colors
      const chainColors = new ColorGenerator();
      chains.forEach((residues, chainID) => {
        const chainGroup = new THREE.Group();
        const chainColor = chainColors.getNextColor();

        // Create backbone trace for this chain
        const backbonePoints = [];
        let prevCA = null;

        residues.forEach((atoms, resSeq) => {
          // Find alpha carbon (CA) for backbone trace
          const ca = atoms.find(a => a.name === 'CA');
          if (ca) {
            backbonePoints.push(new THREE.Vector3(
              ca.x + centerOffset.x,
              ca.y + centerOffset.y,
              ca.z + centerOffset.z
            ));

            if (prevCA) {
              // Create backbone segment
              const start = new THREE.Vector3(
                prevCA.x + centerOffset.x,
                prevCA.y + centerOffset.y,
                prevCA.z + centerOffset.z
              );
              const end = new THREE.Vector3(
                ca.x + centerOffset.x,
                ca.y + centerOffset.y,
                ca.z + centerOffset.z
              );

              const backboneGeom = createCylinderBetweenPoints(start, end, 0.1);
              const backboneMaterial = new THREE.MeshPhongMaterial({
                color: chainColor,
                opacity: 0.7,
                transparent: true
              });
              const backbone = new THREE.Mesh(backboneGeom, backboneMaterial);
              chainGroup.add(backbone);
            }
            prevCA = ca;
          }

          // Create residue group
          const residueGroup = new THREE.Group();

          // Add atoms for this residue
          atoms.forEach(atom => {
            const atomMesh = createAtomMesh(atom, centerOffset);
            residueGroup.add(atomMesh);

            // Create bonds within residue
            atoms.forEach(otherAtom => {
              if (otherAtom.serial > atom.serial) {
                const bondLength = calculateDistance(atom, otherAtom);
                if (bondLength < 2.0) { // Typical bond length threshold
                  const bond = createBond(atom, otherAtom, centerOffset);
                  if (bond) residueGroup.add(bond);
                }
              }
            });
          });

          chainGroup.add(residueGroup);
        });

        // Create smooth backbone curve
        if (backbonePoints.length > 2) {
          const curve = new THREE.CatmullRomCurve3(backbonePoints);
          const tubeGeometry = new THREE.TubeGeometry(curve, backbonePoints.length * 4, 0.2, 8, false);
          const tubeMaterial = new THREE.MeshPhongMaterial({
            color: chainColor,
            opacity: 0.5,
            transparent: true
          });
          const tubeMesh = new THREE.Mesh(tubeGeometry, tubeMaterial);
          chainGroup.add(tubeMesh);
        }

        structureGroup.add(chainGroup);
      });

      // Scale the entire structure
      const size = Math.max(maxX - minX, maxY - minY, maxZ - minZ);
      const scale = 40 / size; // Adjust this value to change overall size
      structureGroup.scale.setScalar(scale);


      const boundingBox = new THREE.Box3().setFromObject(structureGroup);
      const actualCenter = new THREE.Vector3();
      boundingBox.getCenter(actualCenter);
      console.log('Structure bounds:', {
        min: boundingBox.min,
        max: boundingBox.max,
        center: boundingBox.getCenter(new THREE.Vector3())
      });

      console.log('Actual center after creation:', actualCenter);

      // Approach 1: Move structure to exact center
      structureGroup.position.sub(actualCenter);

      // Approach 2: Adjust camera and controls to focus on actual center
      // Calculate proper camera distance based on bounding box size
      const boundingSize = new THREE.Vector3();
      boundingBox.getSize(boundingSize);
      const maxDim = Math.max(boundingSize.x, boundingSize.y, boundingSize.z);
      const fov = camera.fov * (Math.PI / 180);
      const cameraDistance = Math.abs(maxDim / Math.sin(fov / 2));

      // Position camera to look at actual center
      camera.position.set(0, 0, cameraDistance);
      camera.lookAt(actualCenter);

      scene.add(structureGroup);
      return {
        structureGroup,
        center: actualCenter,
        boundingBox,
        updateCamera: (controls) => {
          camera.lookAt(actualCenter);
          if (controls) {
            controls.target.copy(actualCenter);
            controls.update();
          }
        }
      };
    };

    const init = async () => {
      // Scene setup
      scene = new THREE.Scene();
      scene.background = new THREE.Color(0x777777);

      camera = new THREE.PerspectiveCamera(
        75,
        window.innerWidth / window.innerHeight,
        0.1,
        1000
      );
      camera.position.set(0, 0, 100);
      camera.lookAt(0, 0, 0);

      let lastRenderTime = 0;
      const targetFPS = 30; // lower FPS for better performance
      const frameInterval = 1000 / targetFPS;

      renderer = new THREE.WebGLRenderer({ antialias: true, powerPreference: 'high-performance' });
      renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2)); // limit pixel ratio
      renderer.setSize(800, 600);
      mountRef.current.appendChild(renderer.domElement);

      // Lights
      const ambientLight = new THREE.AmbientLight(0x404040);
      scene.add(ambientLight);

      const directionalLight = new THREE.DirectionalLight(0xffffff, 1);
      directionalLight.position.set(10, 10, 10);
      scene.add(directionalLight);

      // Fetch and create genome visualization
      const [data, pdbData] = await Promise.all([fetchGenomeData(), fetchPDBData('7YX4')]);
      if (data && viewMode === 'genome') {
        const { genomeGroup } = createGenomeVisualization(data);
        genomeObject = genomeGroup;

        genomeObject.visible = viewMode === 'genome' || viewMode === 'combined';
      }

      const controls = setupControls();
      if (pdbData && viewMode === 'structure') {
        const structure = parsePDBStructure(pdbData);
        const { structureGroup, updateCamera, center } = createPDBVisualization(structure, scene);
        structureObject = structureGroup
        structureObject.visible = viewMode === 'structure' || viewMode === 'combined';

        // Position structure relative to genome
        structureObject.position.set(50, 0, 0); // Adjust position as needed
      }

      setLoading(false);

      //const animate = () => {
      //  requestAnimationFrame(animate);
      //
      //  controls.update();
      //
      //  if (genomeObject && !controls.enabled) {
      //    genomeObject.rotation.z += 0.001;
      //  }
      //
      //  renderer.render(scene, camera);
      //};
      //
      const animate = (currentTime) => {
        if (currentTime - lastRenderTime > frameInterval) {
          if (controls) controls.update();
          //visualizationInfo.update();
          renderer.render(scene, camera);
          lastRenderTime = currentTime;
        }
        requestAnimationFrame(animate);
      };

      animate();
    };

    init();

    //const handleResize = () => {
    //  const width = mountRef.current.clientWidth;
    //  const height = mountRef.current.clientHeight;
    //
    //  camera.aspect = width / height;
    //  camera.updateProjectionMatrix();
    //  renderer.setSize(width, height);
    //  controls.update();
    //};
    //
    //window.addEventListener('resize', handleResize);

    return () => {
      //window.removeEventListener('resize', handleResize);
      //controls.dispose();
      if (renderer) {
        mountRef.current.removeChild(renderer.domElement);
      }
    };
  }, [viewMode]);

  return (
    <div className="relative h-full m-auto">
      <div ref={mountRef} className="flex items-center justify-center" />

      {loading ? (
        <div className="absolute top-4 left-4 bg-white bg-opacity-75 p-4 rounded shadow">
          Loading genome data...
        </div>
      ) : (
        <>
          <div className="absolute top-4 left-4 bg-white bg-opacity-75 p-4 rounded shadow">
            <h2 className="text-lg font-bold">Mimivirus Genome</h2>
            <p className="text-sm">Green: Genes</p>
            <p className="text-sm">Red: Coding sequences</p>
            <p className="text-sm">Blue: Other features</p>
          </div>

          <h2 className="text-lg font-bold mb-2">View Controls</h2>
          <div className="flex gap-2">
            <button
              className={`px-3 py-1 rounded ${viewMode === 'genome' ? 'bg-blue-500 text-white' : 'bg-gray-200'}`}
              onClick={() => setViewMode('genome')}
            >
              Genome
            </button>
            <button
              className={`px-3 py-1 rounded ${viewMode === 'structure' ? 'bg-blue-500 text-white' : 'bg-gray-200'}`}
              onClick={() => setViewMode('structure')}
            >
              Structure
            </button>
          </div>

          <div className="mt-4">
            <h3 className="font-bold">Legend</h3>
            <p className="text-sm">Green: Genes</p>
            <p className="text-sm">Red: Coding sequences</p>
            <p className="text-sm">Blue: RNA features</p>
          </div>

          {selectedFeature && (
            <div className="absolute bottom-4 left-4 bg-white bg-opacity-75 p-4 rounded shadow max-w-md">
              <h3 className="font-bold">{selectedFeature.type}</h3>
              <p className="text-sm">Location: {selectedFeature.location}</p>
              {selectedFeature.properties.product && (
                <p className="text-sm">
                  Product: {selectedFeature.properties.product}
                </p>
              )}
              {selectedFeature.properties.note && (
                <p className="text-sm">
                  Note: {selectedFeature.properties.note}
                </p>
              )}
            </div>
          )}
        </>
      )
      }
    </div >
  );
};

export default GenomeViewer;
