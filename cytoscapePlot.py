import dash
import dash_cytoscape as cyto
from dash.dependencies import Input, Output
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
cyto.load_extra_layouts()
app = dash.Dash(__name__)
server = app.server

def add_individual_protein(df, source, elements):
    highest = df["Difference"].max()
    n = 0
    for i, r in df.iterrows():
        if n < 15:
            opacity = r["Difference"]/highest
            elements.append({'data': {'id': r["Gene.names"], 'label': r["Gene.names"], 'color': f"rgba(136, 86, 167,{opacity})", "opacity": opacity}, 'classes': 'protein'})
            elements.append(
                {'data': {'source': source, 'target': r["Gene.names"], 'color': f"rgba(136, 86, 167,{opacity})", "opacity": opacity}, 'classes': 'protein-edge'},)
        else:
            break
        n += 1


def add_groups_enriched(edf, elements):
    edf = edf.sort_values(by="Difference", ascending=False)
    golgi = edf[(edf["Golgi"] == "+")]
    #golgi = edf[(edf["C: Golgi"] == "+")]
    golgi_count = len(golgi.index)
    print(golgi_count)
    glyco = golgi[golgi["Glycosylation"] == "+"]
    #glyco = golgi[golgi["Glycosylation genes"] == "+"]
    glyco_count = len(glyco.index)
    print(glyco_count)
    phospha = golgi[golgi["Phosphatases"] == "+"]
    phospha_count = len(phospha.index)
    kinases = golgi[(golgi["Kinases"] == "+") | (golgi["Dark.kinase"] == "+")]
    #kinases = golgi[(golgi["Kinases"] == "+") | (golgi["Dark Kinases"] == "+")]
    kinases_count = len(kinases.index)
    ubi = golgi[golgi["Ub.Pathway"] == "+"]
    ubi_count = len(ubi.index)

    l = [

        {'data': {'id': 'enriched-golgi', 'label': f'Golgi: {golgi_count}', "size": golgi_count * block},
         'classes': 'golgi enriched'},
        #{'data': {'source': 'significant', 'target': 'enriched-golgi'}, 'classes': 'significant-edge'},
        {'data': {'id': 'enriched-glyco', 'label': f'Glycosylation genes: {glyco_count}',
                  "size": glyco_count * block}, 'classes': 'golgi enriched'},
        {'data': {'id': 'enriched-phospha', 'label': f'Phosphatases: {phospha_count}',
                  "size": phospha_count * block}, 'classes': 'golgi enriched'},
        {'data': {'id': 'enriched-kinase', 'label': f'Kinases: {kinases_count}', "size": kinases_count * block},
         'classes': 'golgi enriched'},
        {'data': {'id': 'ubi', 'label': f'Ubiquitin components: {ubi_count}', "size": ubi_count * block},
         'classes': 'golgi enriched'},
        {'data': {'source': 'enriched-golgi', 'target': 'enriched-glyco'}, 'classes': 'golgi-edge enriched'},
        {'data': {'source': 'enriched-golgi', 'target': 'enriched-phospha'},
         'classes': 'golgi-edge enriched'},
        {'data': {'source': 'enriched-golgi', 'target': 'enriched-kinase'},
         'classes': 'golgi-edge enriched'},
        {'data': {'source': 'enriched-golgi', 'target': 'ubi'},
         'classes': 'golgi-edge enriched'},

    ]
    for i in l:
        elements.append(i)
    add_individual_protein(glyco, "enriched-glyco", elements)
    add_individual_protein(phospha, "enriched-phospha", elements)
    add_individual_protein(kinases, "enriched-kinase", elements)
    add_individual_protein(ubi, "ubi", elements)

def add_groups_not_enriched(edf, elements):
    edf = edf.sort_values(by="Difference", ascending=False)
    golgi = edf[(edf["Golgi"] != "+")]
    #golgi = edf[(edf["C: Golgi"] != "+")]
    golgi_count = len(golgi.index)
    print(golgi_count)
    glyco = golgi[golgi["Glycosylation"] == "+"]
    #glyco = golgi[golgi["Glycosylation genes"] == "+"]
    glyco_count = len(glyco.index)
    print(glyco_count)
    phospha = golgi[golgi["Phosphatases"] == "+"]
    phospha_count = len(phospha.index)
    kinases = golgi[(golgi["Kinases"] == "+") | (golgi["Dark.kinase"] == "+")]
    #kinases = golgi[(golgi["Kinases"] == "+") | (golgi["Dark Kinases"] == "+")]
    kinases_count = len(kinases.index)
    ubi = golgi[golgi["Ub.Pathway"] == "+"]
    ubi_count = len(ubi.index)

    l = [

        {'data': {'id': 'non-enriched-golgi', 'label': f'Non-golgi: {golgi_count}', "size": golgi_count * block}, 'classes': 'not-golgi not-enriched'},
        #{'data': {'source': 'significant', 'target': 'non-enriched-golgi'}, 'classes': 'not-golgi significant-edge'},
        {'data': {'id': 'non-enriched-glyco', 'label': f'Glycosylation genes: {glyco_count}',
                  "size": glyco_count * block}, 'classes': 'not-golgi not-enriched'},
        {'data': {'id': 'non-enriched-phospha', 'label': f'Phosphatases: {phospha_count}',
                  "size": phospha_count * block}, 'classes': 'not-golgi not-enriched'},
        {'data': {'id': 'non-enriched-kinase', 'label': f'Kinases: {kinases_count}', "size": kinases_count * block},
         'classes': 'not-golgi not-enriched'},
        {'data': {'id': 'non-enriched-ubi', 'label': f'Ubiquitin components: {ubi_count}', "size": ubi_count * block},
         'classes': 'not-golgi not-enriched'},
        {'data': {'source': 'non-enriched-golgi', 'target': 'non-enriched-glyco'}, 'classes': 'not-golgi-edge not-enriched'},
        {'data': {'source': 'non-enriched-golgi', 'target': 'non-enriched-phospha'}, 'classes': 'not-golgi-edge not-enriched'},
        {'data': {'source': 'non-enriched-golgi', 'target': 'non-enriched-kinase'}, 'classes': 'not-golgi-edge not-enriched'},
        {'data': {'source': 'non-enriched-golgi', 'target': 'non-enriched-ubi'},
         'classes': 'not-golgi-edge not-enriched'},


    ]
    for i in l:
        elements.append(i)
    add_individual_protein(glyco, "non-enriched-glyco", elements)
    add_individual_protein(phospha, "non-enriched-phospha", elements)
    add_individual_protein(kinases, "non-enriched-kinase", elements)
    add_individual_protein(ubi, "non-enriched-ubi", elements)

block = 0.2

#df = pd.read_csv(r"C:\Users\toanp\Downloads\All enriched_For Network.txt", sep="\t")
#df = pd.read_csv(r"C:\Users\toanp\Downloads\GT-IP_Mock-IP_tTest.txt", sep="\t")
df = pd.read_csv(r"C:\Users\toanp\Downloads\GT-IP_WCL_tTest.txt", sep="\t")

df = df[(df["Significant"]=="+")&(df["Difference"] >= 1)]

elements = [
    #{'data': {'id': 'significant', 'label': f'Significant: {len(df.index)}', "size": len(df.index) * block}, 'classes': 'significant'},
]

add_groups_enriched(df, elements)
add_groups_not_enriched(df, elements)

app.layout = html.Div([
    cyto.Cytoscape(
        id='cytoscape',
        elements=elements,
        layout={'name': 'cose', 'idealEdgeLength': 20},
        style={'width': '2000px', 'height': '2000px'},
        stylesheet=[
            {
                'selector': '.significant',
                'style': {
                    'shape': 'ellipse',
                    'background-color': 'rgb(173, 218, 226)',
                }
            },
            {
                'selector': '.not-golgi',
                'style': {
                    'shape': 'ellipse',
                    'background-color': 'rgb(255, 154, 162)',
                }
            },
            {
                'selector': '.not-golgi-edge',
                'style': {
                    'curve-style': 'straight-triangle',
                    "width": 5,
                    'line-color': 'rgb(255, 154, 162)',
                }
            },
            {
                'selector': '.golgi',
                'style': {
                    'shape': 'ellipse',
                    'background-color': 'rgb(255, 218, 193)',
                }
            },
            {
                'selector': '.golgi-edge',
                'style': {
                    'curve-style': 'straight-triangle',
                    "width": 5,
                    'line-color': 'rgb(255, 218, 193)',
                }
            },
            {
                'selector': '.protein',
                'style': {
                    'shape': 'ellipse',
                    'background-color': 'data(color)',
                    'background-opacity': 'data(opacity)',
                    'line-color': 'black'
                }
            },
            {
                'selector': '.protein-edge',
                'style': {
                    'line-color': 'data(color)',
                    'opacity': 'data(opacity)',
                }
            },
            {
                'selector': 'node',
                'style': {
                    "content": "data(label)",
                    "width": "data(size)",
                    "height": "data(size)",
                }
            },
            {
                'selector': '.enriched',
                'style': {
                    'shape': 'ellipse',
                    'background-color': 'rgb(255, 154, 162)',
                    'line-color': 'rgb(255, 154, 162)',
                }
            },
            {
                'selector': '.not-enriched',
                'style': {
                    'shape': 'ellipse',
                    'background-color': 'rgb(255, 218, 193)',
                    'line-color': 'rgb(255, 218, 193)',
                }
            }
        ]
    ),
    html.Div([html.Button("as svg", id="btn-get-svg")])
])
print(elements)
@app.callback(
    Output('image-text', 'children'),
    Input('cytoscape', 'imageData'),
    )
def put_image_string(data):
    return data


@app.callback(
    Output("cytoscape", "generateImage"),
    [
        Input("btn-get-svg", "n_clicks"),
    ])
def get_image(get_svg_clicks):

    # File type to output of 'svg, 'png', 'jpg', or 'jpeg' (alias of 'jpg')


    # 'store': Stores the image data in 'imageDataf' !only jpg/png are supported
    # 'download'`: Downloads the image as a file with all data handling
    # 'both'`: Stores image data and downloads image as file.


    ctx = dash.callback_context
    if ctx.triggered:
        input_id = ctx.triggered[0]["prop_id"].split(".")[0]

        if input_id != "tabs":
            action = "download"
            ftype = input_id.split("-")[-1]
            return {
                'type': 'svg',
                'action': 'download'
            }

    return {
        'type': 'png',
        'action': 'store'
        }

if __name__ == "__main__":
    app.run_server(debug=True)