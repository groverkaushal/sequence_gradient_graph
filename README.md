# Sequence_gradient_graph

This project features a flask web app.

The website takes multiple DNA sequences and some parameters as input.
It then aligns the sequences by multiple sequence alignment using MUSCLE API. Then processess the data according to the parameters given.

The plots are displayed as output.


## Installation

Clone the repsitory
```
git clone https://github.com/Blassreiter0/sequence_gradient_graph.git
```

Create a virtual environment:
```
python -m venv venv
```

Activate the virtual environment:
```
source venv/bin/activate
```

Install the requirements:
```
pip install -r requirements.txt
```


## Usage

Start the Flask app:
```
python rest_app.py
```

Open a web browser and go to http://localhost:5000.
You should see the home page of the Flask app.
